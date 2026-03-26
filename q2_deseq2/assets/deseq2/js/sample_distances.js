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

  const parseRowLabel = (label) => {
    const normalized = String(label || "");
    const parts = normalized
      .split("|")
      .map((part) => part.trim())
      .filter((part) => part !== "");

    const sampleId = parts[0] || normalized;
    const metadata = parts.slice(1).map((part) => {
      const separatorIndex = part.indexOf("=");
      if (separatorIndex === -1) {
        return {
          field: "",
          value: part,
          raw: part,
        };
      }
      return {
        field: part.slice(0, separatorIndex).trim(),
        value: part.slice(separatorIndex + 1).trim(),
        raw: part,
      };
    });

    return {
      raw: normalized,
      sampleId,
      metadata,
    };
  };

  const buildMetadataColorMap = (parsedRowLabels) => {
    const valuesByField = new Map();
    parsedRowLabels.forEach(({ metadata }) => {
      metadata.forEach(({ field, value }) => {
        const fieldKey = field || "_metadata";
        if (!valuesByField.has(fieldKey)) {
          valuesByField.set(fieldKey, []);
        }
        if (!valuesByField.get(fieldKey).includes(value)) {
          valuesByField.get(fieldKey).push(value);
        }
      });
    });

    const colorsByField = new Map();
    Array.from(valuesByField.entries()).forEach(([field, values], fieldIndex) => {
      const colorsByValue = new Map();
      values.forEach((value, valueIndex) => {
        const hue = (fieldIndex * 97 + valueIndex * 41) % 360;
        const saturation = 72;
        const lightness = 88 - (valueIndex % 3) * 7;
        colorsByValue.set(value, {
          background: `hsl(${hue} ${saturation}% ${lightness}%)`,
          border: `hsl(${hue} ${Math.max(42, saturation - 18)}% ${Math.max(52, lightness - 18)}%)`,
          text: `hsl(${hue} 48% 22%)`,
        });
      });
      colorsByField.set(field, colorsByValue);
    });

    return colorsByField;
  };

  const renderRowLabels = (
    parsedRowLabels,
    metadataColors,
    plotHeight,
    rowHeight
  ) => {
    rowLabelNode.innerHTML = "";
    rowLabelNode.style.gridTemplateRows = `repeat(${parsedRowLabels.length}, ${rowHeight}px)`;
    rowLabelNode.style.paddingTop = `${Math.max(0, (plotHeight - parsedRowLabels.length * rowHeight) / 2)}px`;

    parsedRowLabels.forEach((row) => {
      const item = document.createElement("div");
      item.className = "sample-distance-row-label";
      item.style.height = `${rowHeight}px`;
      item.title = row.raw;

      const sample = document.createElement("span");
      sample.className = "sample-distance-row-sample";
      sample.textContent = row.sampleId;
      item.append(sample);

      row.metadata.forEach((entry) => {
        const separator = document.createElement("span");
        separator.className = "sample-distance-row-separator";
        separator.textContent = "|";
        item.append(separator);

        const badge = document.createElement("span");
        badge.className = "sample-distance-row-meta";
        badge.textContent = entry.value || entry.raw;
        badge.title = entry.raw;

        const fieldKey = entry.field || "_metadata";
        const color = metadataColors.get(fieldKey)?.get(entry.value);
        if (color) {
          badge.style.setProperty("--sample-distance-meta-bg", color.background);
          badge.style.setProperty("--sample-distance-meta-border", color.border);
          badge.style.setProperty("--sample-distance-meta-text", color.text);
        }
        item.append(badge);
      });

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
    const parsedRowLabels = rowLabels.map(parseRowLabel);
    const metadataColors = buildMetadataColorMap(parsedRowLabels);
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
      renderRowLabels(
        parsedRowLabels,
        metadataColors,
        plotSize.height,
        plotSize.rowHeight
      );

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
