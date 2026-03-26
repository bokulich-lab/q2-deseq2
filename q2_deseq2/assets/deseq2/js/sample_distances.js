document.addEventListener("DOMContentLoaded", async () => {
  const container = document.getElementById("sample-distance-view");
  const layout = container ? container.parentElement : null;
  const rowLabelNode = document.getElementById("sample-distance-row-labels");
  const legendNode = document.getElementById("sample-distance-legend");
  const specNode = document.getElementById("vega_sample_distance_spec");
  const dataPathNode = document.getElementById("sample_distances_data_path");
  const orderNode = document.getElementById("sample_distance_order");
  const errorNode = document.getElementById("sample-distance-error");

  if (
    !container || !layout || !rowLabelNode || !legendNode
    || !specNode || !dataPathNode || !orderNode
  ) {
    return;
  }

  const parseJsonNode = (node, fallback) => {
    try {
      return JSON.parse(node.textContent);
    } catch {
      return fallback;
    }
  };

  const measureWidth = (element) => {
    return Math.floor(element.getBoundingClientRect().width || 0);
  };

  const clamp = (value, minimum, maximum) => {
    return Math.min(Math.max(value, minimum), maximum);
  };

  const formatLegendValue = (value) => {
    if (!Number.isFinite(value)) {
      return "NA";
    }
    if (Math.abs(value) >= 10) {
      return value.toFixed(1).replace(/\.0$/, "");
    }
    return value.toFixed(2).replace(/0+$/, "").replace(/\.$/, "");
  };

  const computePlotSize = (availableWidth, sampleCount) => {
    const plotWidth = Math.max(240, availableWidth - 8);
    const cellWidth = plotWidth / Math.max(sampleCount, 1);
    const rowHeight = clamp(Math.round(cellWidth * 0.72), 28, 52);
    return {
      width: plotWidth,
      height: Math.max(240, sampleCount * rowHeight),
      rowHeight,
    };
  };

  const renderLegend = (cells) => {
    const distances = cells
      .map((cell) => Number(cell.distance))
      .filter((value) => Number.isFinite(value));
    if (distances.length === 0) {
      legendNode.hidden = true;
      legendNode.innerHTML = "";
      return;
    }

    const maximum = Math.max(...distances);
    const midpoint = maximum / 2;
    legendNode.hidden = false;
    legendNode.innerHTML = `
      <div class="sample-distance-legend-title">Sample distance</div>
      <div class="sample-distance-legend-bar"></div>
      <div class="sample-distance-legend-scale">
        <span>${formatLegendValue(0)}</span>
        <span>${formatLegendValue(midpoint)}</span>
        <span>${formatLegendValue(maximum)}</span>
      </div>
    `;
  };

  const renderRowLabels = (rowLabels, plotHeight, rowHeight) => {
    rowLabelNode.innerHTML = "";
    rowLabelNode.style.gridTemplateRows = `repeat(${rowLabels.length}, ${rowHeight}px)`;
    rowLabelNode.style.paddingTop = `${Math.max(0, (plotHeight - rowLabels.length * rowHeight) / 2)}px`;

    rowLabels.forEach((label) => {
      const item = document.createElement("div");
      item.className = "sample-distance-row-label";
      item.style.height = `${rowHeight}px`;
      item.textContent = label;
      item.title = label;
      rowLabelNode.append(item);
    });
  };

  const showError = (error) => {
    if (errorNode) {
      errorNode.hidden = false;
      errorNode.textContent =
        "The sample-distance heatmap could not be loaded.\n\n" + String(error);
    }
  };

  try {
    const baseSpec = parseJsonNode(specNode, null);
    const sampleOrder = parseJsonNode(orderNode, []);
    const dataPath = parseJsonNode(dataPathNode, "");
    if (
      !baseSpec || !Array.isArray(sampleOrder) || sampleOrder.length === 0 || !dataPath
    ) {
      return;
    }

    const response = await fetch(dataPath);
    if (!response.ok) {
      throw new Error(`Failed to load ${dataPath}: ${response.status}`);
    }

    const cells = await response.json();
    const rowLabelBySample = new Map();
    cells.forEach((cell) => {
      if (!rowLabelBySample.has(cell.sample_y)) {
        rowLabelBySample.set(cell.sample_y, cell.sample_y_label || cell.sample_y);
      }
    });

    const rowLabels = sampleOrder.map(
      (sampleId) => rowLabelBySample.get(sampleId) || sampleId
    );
    const hasTooltipMetadata = cells.some((cell) => (
      (typeof cell.sample_x_metadata === "string" && cell.sample_x_metadata !== "")
      || (typeof cell.sample_y_metadata === "string" && cell.sample_y_metadata !== "")
    ));

    renderLegend(cells);

    let embedResult = null;
    let currentWidth = 0;
    let currentHeight = 0;
    let scheduled = false;

    const applySpecConfig = (renderSpec, plotSize) => {
      renderSpec.width = plotSize.width;
      renderSpec.height = plotSize.height;
      renderSpec.padding = {
        left: 14,
        right: 18,
        top: 8,
        bottom: 84,
      };
      renderSpec.data[0].values = cells;
      renderSpec.legends = [];
      renderSpec.axes = renderSpec.axes.filter((axis) => axis.scale === "x");

      const xScale = renderSpec.scales.find((scale) => scale.name === "x");
      const yScale = renderSpec.scales.find((scale) => scale.name === "y");
      if (xScale) {
        xScale.domain = sampleOrder;
        xScale.paddingInner = 0;
        xScale.paddingOuter = 0;
      }
      if (yScale) {
        yScale.domain = sampleOrder;
        yScale.paddingInner = 0;
        yScale.paddingOuter = 0;
      }

      const xAxis = renderSpec.axes.find((axis) => axis.scale === "x");
      if (xAxis) {
        xAxis.labelFontSize = sampleOrder.length >= 20 ? 10 : 11;
        xAxis.labelAngle = 90;
        xAxis.labelAlign = "left";
        xAxis.labelBaseline = "middle";
        xAxis.title = null;
      }

      if (renderSpec.marks?.[0]?.encode?.update?.y) {
        renderSpec.marks[0].encode.update.y.field = "sample_y";
      }

      if (renderSpec.marks?.[0]?.encode?.enter) {
        renderSpec.marks[0].encode.enter.tooltip = hasTooltipMetadata
          ? {
              signal: "{'Sample A': datum.sample_x, 'Sample A metadata': datum.sample_x_metadata || 'NA', 'Sample B': datum.sample_y, 'Sample B metadata': datum.sample_y_metadata || 'NA', 'Distance': datum.distance == null ? 'NA' : format(datum.distance, '.4f')}"
            }
          : {
              signal: "{'Sample A': datum.sample_x, 'Sample B': datum.sample_y, 'Distance': datum.distance == null ? 'NA' : format(datum.distance, '.4f')}"
            };
      }
    };

    const renderOrResize = async () => {
      const availableWidth = measureWidth(container);
      if (availableWidth < 240) {
        return;
      }

      const plotSize = computePlotSize(availableWidth, sampleOrder.length);
      renderRowLabels(rowLabels, plotSize.height, plotSize.rowHeight);

      if (
        embedResult
        && plotSize.width === currentWidth
        && plotSize.height === currentHeight
      ) {
        return;
      }

      currentWidth = plotSize.width;
      currentHeight = plotSize.height;

      if (!embedResult) {
        const renderSpec = JSON.parse(JSON.stringify(baseSpec));
        applySpecConfig(renderSpec, plotSize);
        embedResult = await vegaEmbed("#sample-distance-view", renderSpec, {
          renderer: "canvas",
          actions: {
            export: { png: true, svg: true },
            source: false,
            compiled: false,
            editor: false,
          },
          downloadFileName: "deseq2-sample-distance-heatmap",
          scaleFactor: { png: 2, svg: 1 },
        });
        return;
      }

      embedResult.view.width(currentWidth);
      embedResult.view.height(currentHeight);
      await embedResult.view.runAsync();
    };

    const scheduleResize = () => {
      if (scheduled) {
        return;
      }
      scheduled = true;
      window.requestAnimationFrame(async () => {
        scheduled = false;
        try {
          await renderOrResize();
        } catch (error) {
          showError(error);
        }
      });
    };

    scheduleResize();
    window.addEventListener("resize", scheduleResize);
    if (typeof ResizeObserver === "function") {
      const observer = new ResizeObserver(() => {
        scheduleResize();
      });
      observer.observe(layout);
      observer.observe(container);
    }
  } catch (error) {
    showError(error);
  }
});
