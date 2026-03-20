document.addEventListener("DOMContentLoaded", async () => {
  const volcanoContainer = document.getElementById("volcano-view");
  const maContainer = document.getElementById("ma-view");
  const volcanoLayout = volcanoContainer ? volcanoContainer.parentElement : null;
  const maLayout = maContainer ? maContainer.parentElement : null;
  const volcanoSpecNode = document.getElementById("vega_volcano_spec");
  const maSpecNode = document.getElementById("vega_ma_spec");
  const volcanoDataPathNode = document.getElementById("volcano_data_path");
  const maDataPathNode = document.getElementById("ma_data_path");
  const volcanoErrorNode = document.getElementById("volcano-error");
  const maErrorNode = document.getElementById("ma-error");
  const volcanoFallbackNode = document.getElementById("volcano-fallback");
  const maFallbackNode = document.getElementById("ma-fallback");
  const alphaInput = document.getElementById("alpha-cutoff-input");
  const lfcInput = document.getElementById("lfc-cutoff-input");
  const resetButton = document.getElementById("plot-controls-reset");
  const significantNode = document.getElementById("summary-significant-features");

  if (
    !volcanoContainer || !maContainer || !volcanoLayout || !maLayout
    || !volcanoSpecNode || !maSpecNode
    || !volcanoDataPathNode || !maDataPathNode
    || !alphaInput || !lfcInput
  ) {
    return;
  }

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

  const alphaBounds = readBounds(alphaInput, 0.05);
  const lfcBounds = readBounds(lfcInput, 1.0);
  let currentAlpha = clamp(alphaBounds.defaultValue, alphaBounds.min, alphaBounds.max);
  let currentLfc = clamp(lfcBounds.defaultValue, lfcBounds.min, lfcBounds.max);

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
    const [volcanoResponse, maResponse] = await Promise.all([
      fetch(volcanoDataPathNode.textContent.trim()),
      fetch(maDataPathNode.textContent.trim()),
    ]);
    if (!volcanoResponse.ok) {
      throw new Error(`Failed to load volcano data: ${volcanoResponse.status}`);
    }
    if (!maResponse.ok) {
      throw new Error(`Failed to load MA data: ${maResponse.status}`);
    }

    const [volcanoData, maData] = await Promise.all([
      volcanoResponse.json(),
      maResponse.json(),
    ]);

    const volcanoSpec = JSON.parse(volcanoSpecNode.textContent);
    volcanoSpec.data[0].values = volcanoData;
    volcanoSpec.width = measureWidth(volcanoLayout);

    const maSpec = JSON.parse(maSpecNode.textContent);
    maSpec.data[0].values = maData;
    maSpec.width = measureWidth(maLayout);

    const [volcanoResult, maResult] = await Promise.all([
      vegaEmbed("#volcano-view", volcanoSpec, { actions: false, renderer: "canvas" }),
      vegaEmbed("#ma-view", maSpec, { actions: false, renderer: "canvas" }),
    ]);

    if (volcanoFallbackNode) {
      volcanoFallbackNode.hidden = true;
    }
    if (maFallbackNode) {
      maFallbackNode.hidden = true;
    }

    const views = [
      { view: volcanoResult.view, layout: volcanoLayout, width: volcanoSpec.width },
      { view: maResult.view, layout: maLayout, width: maSpec.width },
    ];

    const plottableData = volcanoData.filter((datum) => (
      Number.isFinite(datum.log2FoldChange)
      && Number.isFinite(datum.padj)
      && datum.padj > 0
    ));

    const updateSummary = (alphaCutoff, lfcCutoff) => {
      const significantCount = plottableData.filter((datum) => (
        datum.padj <= alphaCutoff && Math.abs(datum.log2FoldChange) >= lfcCutoff
      )).length;
      if (significantNode) {
        significantNode.textContent = significantCount.toLocaleString();
      }
    };

    const pushControlState = () => {
      syncControlInputs();
      updateSummary(currentAlpha, currentLfc);
      views.forEach(({ view }) => {
        view.signal("alphaCutoff", currentAlpha);
        view.signal("lfcCutoff", currentLfc);
      });
      Promise.all(views.map(({ view }) => view.runAsync())).catch(() => {});
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
    pushControlState();

    volcanoResult.view.addSignalListener("hoveredFeatureId", (_, value) => {
      maResult.view.signal("linkedHoveredFeatureId", value).runAsync();
    });
    maResult.view.addSignalListener("hoveredFeatureId", (_, value) => {
      volcanoResult.view.signal("linkedHoveredFeatureId", value).runAsync();
    });

    let resizeFrame = null;
    const resizePlots = () => {
      if (resizeFrame !== null) {
        window.cancelAnimationFrame(resizeFrame);
      }
      resizeFrame = window.requestAnimationFrame(() => {
        resizeFrame = null;
        let needsRun = false;
        views.forEach((plot) => {
          const nextWidth = measureWidth(plot.layout);
          if (nextWidth !== plot.width) {
            plot.width = nextWidth;
            plot.view.width(nextWidth).resize();
            needsRun = true;
          }
        });
        if (needsRun) {
          Promise.all(views.map(({ view }) => view.runAsync())).catch(() => {});
        }
      });
    };

    window.addEventListener("resize", resizePlots);
  } catch (error) {
    if (volcanoErrorNode) {
      volcanoErrorNode.hidden = false;
      volcanoErrorNode.textContent =
        "Interactive volcano plot could not be rendered.\n\n" + String(error);
    }
    if (maErrorNode) {
      maErrorNode.hidden = false;
      maErrorNode.textContent =
        "Interactive MA plot could not be rendered.\n\n" + String(error);
    }
    if (volcanoFallbackNode) {
      volcanoFallbackNode.hidden = false;
    }
    if (maFallbackNode) {
      maFallbackNode.hidden = false;
    }
  }
});
