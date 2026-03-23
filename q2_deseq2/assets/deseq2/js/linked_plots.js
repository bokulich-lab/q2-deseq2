document.addEventListener("DOMContentLoaded", async () => {
  const volcanoContainer = document.getElementById("volcano-view");
  const maContainer = document.getElementById("ma-view");
  const volcanoLayout = volcanoContainer ? volcanoContainer.parentElement : null;
  const maLayout = maContainer ? maContainer.parentElement : null;
  const volcanoSpecNode = document.getElementById("vega_volcano_spec");
  const maSpecNode = document.getElementById("vega_ma_spec");
  const resultsDataPathNode = document.getElementById("results_data_path");
  const effectOptionsNode = document.getElementById("effect_options");
  const defaultEffectIdNode = document.getElementById("default_effect_id");
  const effectSelect = document.getElementById("effect-select");
  const effectLabelNodes = document.querySelectorAll("[data-selected-effect]");
  const volcanoErrorNode = document.getElementById("volcano-error");
  const maErrorNode = document.getElementById("ma-error");
  const volcanoFallbackNode = document.getElementById("volcano-fallback");
  const maFallbackNode = document.getElementById("ma-fallback");
  const alphaInput = document.getElementById("alpha-cutoff-input");
  const lfcInput = document.getElementById("lfc-cutoff-input");
  const resetButton = document.getElementById("plot-controls-reset");
  const significantNode = document.getElementById("summary-significant-features");
  const testedNode = document.getElementById("summary-total-features");
  const annotatedNode = document.getElementById("summary-annotated-features");

  if (
    !volcanoContainer || !maContainer || !volcanoLayout || !maLayout
    || !volcanoSpecNode || !maSpecNode
    || !resultsDataPathNode || !effectOptionsNode || !defaultEffectIdNode
    || !effectSelect || !alphaInput || !lfcInput
  ) {
    return;
  }

  const effectPersistenceKey = "q2_deseq2_selected_effect";
  const legacyComparisonPersistenceKey = "q2_deseq2_selected_comparison";
  const alphaPersistenceKey = "q2_deseq2_alpha_cutoff";
  const lfcPersistenceKey = "q2_deseq2_lfc_cutoff";

  const parseJsonNode = (node, fallback) => {
    try {
      return JSON.parse(node.textContent);
    } catch {
      return fallback;
    }
  };

  const formatComparison = (testLevel, referenceLevel) => {
    if (!testLevel || !referenceLevel) {
      return "";
    }
    return `${testLevel} vs. ${referenceLevel}`;
  };

  const deriveLegacyEffectId = (record, fallbackLabel) => {
    const testLevel = typeof record.test_level === "string" ? record.test_level.trim() : "";
    const referenceLevel = typeof record.reference_level === "string" ? record.reference_level.trim() : "";
    const comparison = typeof record.comparison === "string" ? record.comparison.trim() : "";
    if (testLevel && referenceLevel) {
      return `legacy::${testLevel}::${referenceLevel}`;
    }
    if (comparison) {
      return `legacy::${comparison}`;
    }
    return `legacy::${fallbackLabel || "effect"}`;
  };

  const deriveLegacyEffectLabel = (record, fallbackLabel) => {
    if (typeof record.comparison === "string" && record.comparison.trim() !== "") {
      return record.comparison.trim();
    }
    const comparison = formatComparison(record.test_level, record.reference_level);
    return comparison || fallbackLabel || "Effect";
  };

  const parseNumeric = (value) => {
    if (value === null || value === "") {
      return null;
    }
    if (typeof value === "number") {
      return Number.isFinite(value) ? value : null;
    }
    if (typeof value === "string") {
      const trimmed = value.trim();
      if (trimmed === "") {
        return null;
      }
      const parsed = Number.parseFloat(trimmed);
      return Number.isFinite(parsed) ? parsed : null;
    }
    return null;
  };

  const formatAlpha = (value) => {
    const numeric = Number.parseFloat(value);
    if (!Number.isFinite(numeric)) {
      return "0.05";
    }
    return numeric < 0.001
      ? numeric.toFixed(4)
      : numeric.toFixed(3).replace(/0+$/, "").replace(/\.$/, "");
  };

  const formatLfc = (value) => {
    const numeric = Number.parseFloat(value);
    if (!Number.isFinite(numeric)) {
      return "1.0";
    }
    return numeric.toFixed(1).replace(/\.0$/, "");
  };

  const measureWidth = (element) => {
    const width = Math.floor(element.getBoundingClientRect().width || 0);
    return Math.max(320, width || 860);
  };

  const clamp = (value, minimum, maximum) => {
    return Math.min(Math.max(value, minimum), maximum);
  };

  const readBounds = (input, fallback) => {
    return {
      min: Number.isFinite(Number.parseFloat(input.min))
        ? Number.parseFloat(input.min)
        : fallback,
      max: Number.isFinite(Number.parseFloat(input.max))
        ? Number.parseFloat(input.max)
        : fallback,
      defaultValue: Number.isFinite(Number.parseFloat(input.defaultValue))
        ? Number.parseFloat(input.defaultValue)
        : fallback,
    };
  };

  const readPersistedEffect = (availableEffectIds, effectOptions) => {
    try {
      const persisted = window.localStorage.getItem(effectPersistenceKey);
      if (persisted && availableEffectIds.includes(persisted)) {
        return persisted;
      }
      const legacyPersisted = window.localStorage.getItem(legacyComparisonPersistenceKey);
      if (legacyPersisted) {
        const matchingEffect = effectOptions.find((option) => option.effect_label === legacyPersisted);
        if (matchingEffect) {
          return matchingEffect.effect_id;
        }
      }
    } catch {}
    return null;
  };

  const persistEffect = (effectId) => {
    try {
      window.localStorage.setItem(effectPersistenceKey, effectId);
    } catch {}
  };

  const readPersistedNumeric = (key, bounds) => {
    try {
      const persisted = Number.parseFloat(window.localStorage.getItem(key));
      if (Number.isFinite(persisted)) {
        return clamp(persisted, bounds.min, bounds.max);
      }
    } catch {}
    return null;
  };

  const persistNumeric = (key, value) => {
    try {
      window.localStorage.setItem(key, String(value));
    } catch {}
  };

  const effectOptions = parseJsonNode(effectOptionsNode, []);
  const configuredDefaultEffectId = parseJsonNode(defaultEffectIdNode, "");
  const resultsDataPath = parseJsonNode(resultsDataPathNode, "");

  const alphaBounds = readBounds(alphaInput, 0.05);
  const lfcBounds = readBounds(lfcInput, 1.0);
  let currentAlpha = readPersistedNumeric(alphaPersistenceKey, alphaBounds)
    ?? clamp(alphaBounds.defaultValue, alphaBounds.min, alphaBounds.max);
  let currentLfc = readPersistedNumeric(lfcPersistenceKey, lfcBounds)
    ?? clamp(lfcBounds.defaultValue, lfcBounds.min, lfcBounds.max);
  let activePlotData = [];

  const parseControlValue = (input, bounds, fallback) => {
    const parsed = Number.parseFloat(input.value);
    if (!Number.isFinite(parsed)) {
      return fallback;
    }
    return clamp(parsed, bounds.min, bounds.max);
  };

  const syncControlInputs = () => {
    alphaInput.value = formatAlpha(currentAlpha);
    lfcInput.value = formatLfc(currentLfc);
  };

  try {
    const response = await fetch(resultsDataPath);
    if (!response.ok) {
      throw new Error(`Failed to load ${resultsDataPath}: ${response.status}`);
    }

    const payload = await response.json();
    const columns = payload.columns || [];
    const rows = payload.data || [];
    const records = rows.map((row) => {
      const record = {};
      columns.forEach((name, index) => {
        record[name] = row[index];
      });
      const fallbackLabel = deriveLegacyEffectLabel(record, "Effect");
      record.effect_id = record.effect_id || deriveLegacyEffectId(record, fallbackLabel);
      record.effect_label = record.effect_label || fallbackLabel;
      return record;
    });

    const rowsByEffect = new Map();
    const labelByEffect = new Map();
    effectOptions.forEach((option) => {
      labelByEffect.set(option.effect_id, option.effect_label);
      rowsByEffect.set(option.effect_id, []);
    });
    records.forEach((record) => {
      if (!rowsByEffect.has(record.effect_id)) {
        rowsByEffect.set(record.effect_id, []);
      }
      if (!labelByEffect.has(record.effect_id)) {
        labelByEffect.set(record.effect_id, record.effect_label);
      }
      rowsByEffect.get(record.effect_id).push(record);
    });

    const availableEffectIds = Array.from(rowsByEffect.keys()).filter(Boolean);
    if (availableEffectIds.length === 0) {
      throw new Error("No effects were available in the DESeq2 results table.");
    }

    availableEffectIds.forEach((effectId) => {
      const option = document.createElement("option");
      option.value = effectId;
      option.textContent = labelByEffect.get(effectId) || effectId;
      effectSelect.append(option);
    });

    const initialEffectId = readPersistedEffect(availableEffectIds, effectOptions)
      || (availableEffectIds.includes(configuredDefaultEffectId)
        ? configuredDefaultEffectId
        : availableEffectIds[0]);
    effectSelect.value = initialEffectId;

    const buildPlotData = (effectId) => {
      return (rowsByEffect.get(effectId) || []).map((row) => ({
        feature_id: row.feature_id == null ? "" : String(row.feature_id),
        gene_name: row.gene_name == null || row.gene_name === "" ? null : String(row.gene_name),
        product: row.product == null || row.product === "" ? null : String(row.product),
        log2FoldChange: parseNumeric(row.log2FoldChange),
        pvalue: parseNumeric(row.pvalue),
        padj: parseNumeric(row.padj),
        baseMean: parseNumeric(row.baseMean),
      }));
    };

    const setEffectLabel = (effectId) => {
      const label = labelByEffect.get(effectId) || effectId;
      effectLabelNodes.forEach((node) => {
        node.textContent = label;
      });
    };

    const volcanoSpec = JSON.parse(volcanoSpecNode.textContent);
    volcanoSpec.width = measureWidth(volcanoLayout);
    volcanoSpec.data[0].values = buildPlotData(initialEffectId);

    const maSpec = JSON.parse(maSpecNode.textContent);
    maSpec.width = measureWidth(maLayout);
    maSpec.data[0].values = buildPlotData(initialEffectId);

    const volcanoDataName = Array.isArray(volcanoSpec.data) && volcanoSpec.data[0]
      ? volcanoSpec.data[0].name
      : "points";
    const maDataName = Array.isArray(maSpec.data) && maSpec.data[0]
      ? maSpec.data[0].name
      : "points";

    const [volcanoResult, maResult] = await Promise.all([
      vegaEmbed("#volcano-view", volcanoSpec, {
        renderer: "canvas",
        actions: {
          export: { png: true, svg: true },
          source: false,
          compiled: false,
          editor: false,
        },
        downloadFileName: "deseq2-volcano-plot",
        scaleFactor: { png: 2, svg: 1 },
      }),
      vegaEmbed("#ma-view", maSpec, {
        renderer: "canvas",
        actions: {
          export: { png: true, svg: true },
          source: false,
          compiled: false,
          editor: false,
        },
        downloadFileName: "deseq2-ma-plot",
        scaleFactor: { png: 2, svg: 1 },
      }),
    ]);

    if (volcanoFallbackNode) {
      volcanoFallbackNode.hidden = true;
    }
    if (maFallbackNode) {
      maFallbackNode.hidden = true;
    }

    const views = [
      {
        view: volcanoResult.view,
        layout: volcanoLayout,
        width: volcanoSpec.width,
        dataName: volcanoDataName,
      },
      {
        view: maResult.view,
        layout: maLayout,
        width: maSpec.width,
        dataName: maDataName,
      },
    ];

    const updateSummary = (alphaCutoff, lfcCutoff) => {
      const plottableData = activePlotData.filter((datum) => (
        Number.isFinite(datum.log2FoldChange)
        && Number.isFinite(datum.padj)
        && datum.padj > 0
      ));
      const significantCount = plottableData.filter((datum) => (
        datum.padj <= alphaCutoff && Math.abs(datum.log2FoldChange) >= lfcCutoff
      )).length;
      const annotatedCount = activePlotData.filter((datum) => (
        (typeof datum.gene_name === "string" && datum.gene_name.trim() !== "")
        || (typeof datum.product === "string" && datum.product.trim() !== "")
      )).length;

      if (significantNode) {
        significantNode.textContent = significantCount.toLocaleString();
      }
      if (testedNode) {
        testedNode.textContent = activePlotData.length.toLocaleString();
      }
      if (annotatedNode) {
        annotatedNode.textContent = annotatedCount.toLocaleString();
      }
    };

    const pushControlState = () => {
      syncControlInputs();
      persistNumeric(alphaPersistenceKey, currentAlpha);
      persistNumeric(lfcPersistenceKey, currentLfc);
      updateSummary(currentAlpha, currentLfc);
      views.forEach(({ view }) => {
        view.signal("alphaCutoff", currentAlpha);
        view.signal("lfcCutoff", currentLfc);
      });
      Promise.all(views.map(({ view }) => view.runAsync())).catch(() => {});
    };

    const applyEffect = (effectId) => {
      activePlotData = buildPlotData(effectId);
      views.forEach(({ view, dataName }) => {
        view.change(
          dataName,
          vega
            .changeset()
            .remove(() => true)
            .insert(activePlotData),
        );
      });
      setEffectLabel(effectId);
      effectSelect.value = effectId;
      persistEffect(effectId);
      pushControlState();
    };

    const commitControlState = () => {
      currentAlpha = parseControlValue(alphaInput, alphaBounds, currentAlpha);
      currentLfc = parseControlValue(lfcInput, lfcBounds, currentLfc);
      pushControlState();
    };

    const commitControlInput = (event) => {
      if (!event.target.validity.valid || event.target.value === "") {
        return;
      }
      commitControlState();
    };

    effectSelect.addEventListener("change", () => {
      applyEffect(effectSelect.value);
    });
    alphaInput.addEventListener("input", commitControlInput);
    alphaInput.addEventListener("change", commitControlState);
    alphaInput.addEventListener("blur", commitControlState);
    lfcInput.addEventListener("input", commitControlInput);
    lfcInput.addEventListener("change", commitControlState);
    lfcInput.addEventListener("blur", commitControlState);
    if (resetButton) {
      resetButton.addEventListener("click", () => {
        currentAlpha = clamp(alphaBounds.defaultValue, alphaBounds.min, alphaBounds.max);
        currentLfc = clamp(lfcBounds.defaultValue, lfcBounds.min, lfcBounds.max);
        pushControlState();
      });
    }

    const resizePlots = () => {
      views.forEach((entry) => {
        const width = measureWidth(entry.layout);
        if (width !== entry.width) {
          entry.width = width;
          entry.view.width(width);
          entry.view.runAsync().catch(() => {});
        }
      });
    };
    window.addEventListener("resize", resizePlots);

    applyEffect(initialEffectId);
  } catch (error) {
    if (volcanoErrorNode) {
      volcanoErrorNode.hidden = false;
      volcanoErrorNode.textContent =
        "The volcano plot could not be loaded.\n\n" + String(error);
    }
    if (maErrorNode) {
      maErrorNode.hidden = false;
      maErrorNode.textContent =
        "The MA plot could not be loaded.\n\n" + String(error);
    }
  }
});
