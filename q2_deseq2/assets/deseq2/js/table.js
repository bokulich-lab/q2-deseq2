document.addEventListener("DOMContentLoaded", async () => {
  const dataPathNode = document.getElementById("results_data_path");
  const effectOptionsNode = document.getElementById("effect_options");
  const defaultEffectIdNode = document.getElementById("default_effect_id");
  const effectSelect = document.getElementById("results-effect-select");
  const effectLabelNode = document.getElementById("results-selected-effect");
  const alphaInput = document.getElementById("results-alpha-cutoff-input");
  const lfcInput = document.getElementById("results-lfc-cutoff-input");
  const resetButton = document.getElementById("results-filter-reset");
  const errorNode = document.getElementById("results-table-error");
  const tableNode = document.getElementById("results-table");

  if (
    !dataPathNode || !effectOptionsNode || !defaultEffectIdNode
    || !effectSelect || !effectLabelNode
    || !alphaInput || !lfcInput || !tableNode
  ) {
    return;
  }

  const effectPersistenceKey = "q2_deseq2_selected_effect";
  const legacyComparisonPersistenceKey = "q2_deseq2_selected_comparison";
  const alphaPersistenceKey = "q2_deseq2_alpha_cutoff";
  const lfcPersistenceKey = "q2_deseq2_lfc_cutoff";
  let tableInstance = null;

  const parseJsonNode = (node, fallback) => {
    try {
      return JSON.parse(node.textContent);
    } catch {
      return fallback;
    }
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

  const formatScientific = (value) => {
    if (!Number.isFinite(value)) {
      return "";
    }
    if (value === 0) {
      return "0.00e+0";
    }
    return value.toExponential(2);
  };

  const formatFixed = (value) => {
    if (!Number.isFinite(value)) {
      return "";
    }
    return value.toFixed(2);
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

  const isNumericColumn = (rows, columnIndex) => {
    return rows.every((row) => {
      const value = row[columnIndex];
      return value === null || value === "" || parseNumeric(value) !== null;
    });
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

  const alphaBounds = readBounds(alphaInput, 0.05);
  const lfcBounds = readBounds(lfcInput, 1.0);
  let currentAlpha = readPersistedNumeric(alphaPersistenceKey, alphaBounds)
    ?? clamp(alphaBounds.defaultValue, alphaBounds.min, alphaBounds.max);
  let currentLfc = readPersistedNumeric(lfcPersistenceKey, lfcBounds)
    ?? clamp(lfcBounds.defaultValue, lfcBounds.min, lfcBounds.max);

  const syncControlInputs = () => {
    alphaInput.value = formatAlpha(currentAlpha);
    lfcInput.value = formatLfc(currentLfc);
  };

  const parseControlValue = (input, bounds, fallback) => {
    const parsed = Number.parseFloat(input.value);
    if (!Number.isFinite(parsed)) {
      return fallback;
    }
    return clamp(parsed, bounds.min, bounds.max);
  };

  try {
    const dataPath = parseJsonNode(dataPathNode, "");
    const configuredDefaultEffectId = parseJsonNode(defaultEffectIdNode, "");
    const effectOptions = parseJsonNode(effectOptionsNode, []);

    const response = await fetch(dataPath);
    if (!response.ok) {
      throw new Error(`Failed to load ${dataPath}: ${response.status}`);
    }

    const payload = await response.json();
    const columns = payload.columns || [];
    const rows = payload.data || [];
    const effectIdColumnIndex = columns.indexOf("effect_id");
    const effectLabelColumnIndex = columns.indexOf("effect_label");
    const comparisonColumnIndex = columns.indexOf("comparison");
    const testLevelColumnIndex = columns.indexOf("test_level");
    const referenceLevelColumnIndex = columns.indexOf("reference_level");
    const padjColumnIndex = columns.indexOf("padj");
    const lfcColumnIndex = columns.indexOf("log2FoldChange");
    const rowsByEffect = new Map();
    const labelByEffect = new Map();

    effectOptions.forEach((option) => {
      rowsByEffect.set(option.effect_id, []);
      labelByEffect.set(option.effect_id, option.effect_label);
    });

    rows.forEach((row) => {
      const record = {};
      columns.forEach((name, index) => {
        record[name] = row[index];
      });
      const fallbackLabel = effectLabelColumnIndex >= 0 && row[effectLabelColumnIndex]
        ? row[effectLabelColumnIndex]
        : deriveLegacyEffectLabel(record, "Effect");
      const effectId = effectIdColumnIndex >= 0 && row[effectIdColumnIndex]
        ? row[effectIdColumnIndex]
        : deriveLegacyEffectId(record, fallbackLabel);
      const effectLabel = effectLabelColumnIndex >= 0 && row[effectLabelColumnIndex]
        ? row[effectLabelColumnIndex]
        : fallbackLabel;

      if (!rowsByEffect.has(effectId)) {
        rowsByEffect.set(effectId, []);
      }
      if (!labelByEffect.has(effectId)) {
        labelByEffect.set(effectId, effectLabel);
      }
      rowsByEffect.get(effectId).push(row);
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

    const hiddenColumns = new Set([
      "effect_id",
      "effect_label",
      "effect_kind",
      "effect_expression",
      "comparison",
      "test_level",
      "reference_level",
    ]);
    const columnConfig = columns.map((name, columnIndex) => {
      const numericColumn = isNumericColumn(rows, columnIndex);
      const useFixedDecimal =
        name === "log2FoldChange" || name === "lfcSE" || name === "stat";
      const width =
        name === "feature_id" ? "17.85rem"
          : name === "pvalue" || name === "padj" ? "9rem"
            : name === "gene_name" ? "15rem"
              : null;

      return {
        title: name,
        defaultContent: "",
        visible: !hiddenColumns.has(name),
        width,
        render(data, type) {
          if (!numericColumn || data === null || data === "") {
            return data ?? "";
          }
          const numericValue = parseNumeric(data);
          if (numericValue === null) {
            return data ?? "";
          }
          if (type === "display" || type === "filter") {
            if (useFixedDecimal) {
              return formatFixed(numericValue);
            }
            return formatScientific(numericValue);
          }
          return numericValue;
        },
      };
    });

    let currentEffectId = readPersistedEffect(availableEffectIds, effectOptions)
      || (availableEffectIds.includes(configuredDefaultEffectId)
        ? configuredDefaultEffectId
        : availableEffectIds[0]);

    const rowPassesFilters = (row) => {
      if (padjColumnIndex < 0 || lfcColumnIndex < 0) {
        return true;
      }

      const padjValue = parseNumeric(row[padjColumnIndex]);
      const lfcValue = parseNumeric(row[lfcColumnIndex]);
      return Number.isFinite(padjValue)
        && padjValue > 0
        && padjValue <= currentAlpha
        && Number.isFinite(lfcValue)
        && Math.abs(lfcValue) >= currentLfc;
    };

    const renderRows = (filteredRows) => {
      if (tableInstance === null) {
        tableInstance = new DataTable(tableNode, {
          columns: columnConfig,
          data: filteredRows,
          pageLength: 25,
          scrollX: true,
        });
        return;
      }

      tableInstance.clear();
      tableInstance.rows.add(filteredRows);
      tableInstance.draw(false);
    };

    const applyFilters = () => {
      syncControlInputs();
      persistEffect(currentEffectId);
      persistNumeric(alphaPersistenceKey, currentAlpha);
      persistNumeric(lfcPersistenceKey, currentLfc);

      const effectRows = rowsByEffect.get(currentEffectId) || [];
      const filteredRows = effectRows.filter(rowPassesFilters);
      renderRows(filteredRows);

      effectLabelNode.textContent = labelByEffect.get(currentEffectId) || currentEffectId;
      effectSelect.value = currentEffectId;
    };

    const commitControlState = () => {
      currentAlpha = parseControlValue(alphaInput, alphaBounds, currentAlpha);
      currentLfc = parseControlValue(lfcInput, lfcBounds, currentLfc);
      applyFilters();
    };

    const commitControlInput = (event) => {
      if (!event.target.validity.valid || event.target.value === "") {
        return;
      }
      commitControlState();
    };

    effectSelect.addEventListener("change", () => {
      currentEffectId = effectSelect.value;
      applyFilters();
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
        applyFilters();
      });
    }

    applyFilters();
  } catch (error) {
    if (errorNode) {
      errorNode.hidden = false;
      errorNode.textContent =
        "The result table could not be loaded.\n\n" + String(error);
    }
  }
});
