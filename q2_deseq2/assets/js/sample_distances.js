document.addEventListener("DOMContentLoaded", async () => {
  const parseJsonNode = (node, fallback) => {
    if (!node) {
      return fallback;
    }
    try {
      return JSON.parse(node.textContent);
    } catch {
      return fallback;
    }
  };

  const measureWidth = (element) => {
    if (!element) {
      return 0;
    }
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

  const formatPercent = (value) => {
    if (!Number.isFinite(value)) {
      return "0";
    }
    return value.toFixed(1).replace(/\.0$/, "");
  };

  const heatmapPadding = {
    left: 0,
    right: 8,
    top: 8,
    bottom: 104,
  };

  const countMatrixPadding = {
    left: 0,
    right: 0,
    top: 0,
    bottom: 0,
  };

  const showError = (node, label, error) => {
    if (!node) {
      return;
    }
    node.hidden = false;
    node.textContent = `${label} could not be loaded.\n\n${String(error)}`;
  };

  const scheduleViewResize = (callback, observedElements) => {
    let scheduled = false;

    const schedule = () => {
      if (scheduled) {
        return;
      }
      scheduled = true;
      window.requestAnimationFrame(async () => {
        scheduled = false;
        await callback();
      });
    };

    window.addEventListener("resize", schedule);
    if (typeof ResizeObserver === "function") {
      const observer = new ResizeObserver(() => {
        schedule();
      });
      observedElements.forEach((element) => {
        if (element) {
          observer.observe(element);
        }
      });
    }

    schedule();
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

  const buildGroupPalette = (points, groupField) => {
    if (!groupField) {
      return {
        domain: [],
        colors: new Map(),
      };
    }

    const parsedRows = points.map((point) => ({
      metadata: point.group_value
        ? [{ field: groupField, value: point.group_value }]
        : [],
    }));
    const metadataColors = buildMetadataColorMap(parsedRows);
    const colors = metadataColors.get(groupField) || new Map();
    const domain = [];

    points.forEach((point) => {
      const value = point.group_value || "Samples";
      if (!domain.includes(value)) {
        domain.push(value);
      }
    });

    return {
      domain,
      colors,
    };
  };

  const renderRowLabels = (
    rowLabelNode,
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

  const renderHeatmapLegend = (legendNode, values, title) => {
    const numericValues = values
      .map((value) => Number(value))
      .filter((value) => Number.isFinite(value));
    if (numericValues.length === 0) {
      legendNode.hidden = true;
      legendNode.innerHTML = "";
      return;
    }

    const minimum = Math.min(...numericValues);
    const maximum = Math.max(...numericValues);
    const midpoint = minimum + (maximum - minimum) / 2;
    legendNode.hidden = false;
    legendNode.innerHTML = `
      <div class="sample-distance-legend-title">${title}</div>
      <div class="sample-distance-legend-bar"></div>
      <div class="sample-distance-legend-scale">
        <span>${formatLegendValue(minimum)}</span>
        <span>${formatLegendValue(midpoint)}</span>
        <span>${formatLegendValue(maximum)}</span>
      </div>
    `;
  };

  const renderDistanceLegend = (legendNode, cells) => {
    renderHeatmapLegend(
      legendNode,
      cells.map((cell) => cell.distance),
      "Sample distance"
    );
  };

  const parseMetadataText = (text) => {
    return String(text || "")
      .split(";")
      .map((part) => part.trim())
      .filter((part) => part !== "")
      .map((part) => {
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
  };

  const abbreviateMetadataValue = (value) => {
    const text = String(value || "").trim();
    if (text === "") {
      return "";
    }

    const tokens = text.split(/[^A-Za-z0-9]+/).filter(Boolean);
    if (tokens.length > 1) {
      return tokens
        .map((token) => {
          const head = token.slice(0, 1).toUpperCase();
          const digits = token.replace(/[^0-9]/g, "");
          return `${head}${digits}`;
        })
        .join("")
        .slice(0, 3);
    }

    const alphaNumeric = text.replace(/[^A-Za-z0-9]/g, "");
    if (alphaNumeric.length <= 3) {
      return alphaNumeric.toUpperCase();
    }
    return alphaNumeric.slice(0, 1).toUpperCase();
  };

  const renderCountMatrixSampleLabels = (
    labelNode,
    samples,
    metadataColors,
    plotWidth
  ) => {
    labelNode.innerHTML = "";
    labelNode.style.gridTemplateColumns = `repeat(${samples.length}, minmax(0, 1fr))`;
    labelNode.style.width = `${plotWidth}px`;

    samples.forEach((sample) => {
      const item = document.createElement("div");
      item.className = "count-matrix-sample-header";
      item.title = sample.sample_label || sample.sample_id;

      const sampleName = document.createElement("div");
      sampleName.className = "count-matrix-sample-name";
      sampleName.textContent = sample.sample_id;
      item.append(sampleName);

      const metadataStrip = document.createElement("div");
      metadataStrip.className = "count-matrix-sample-meta-strip";
      parseMetadataText(sample.sample_metadata).forEach((entry) => {
        const badge = document.createElement("span");
        badge.className = "sample-distance-row-meta";
        badge.textContent = abbreviateMetadataValue(entry.value || entry.raw);
        badge.title = entry.raw;

        const fieldKey = entry.field || "_metadata";
        const color = metadataColors.get(fieldKey)?.get(entry.value);
        if (color) {
          badge.style.setProperty("--sample-distance-meta-bg", color.background);
          badge.style.setProperty("--sample-distance-meta-border", color.border);
          badge.style.setProperty("--sample-distance-meta-text", color.text);
        }
        metadataStrip.append(badge);
      });

      item.append(metadataStrip);
      labelNode.append(item);
    });
  };

  const renderSimpleRowLabels = (
    rowLabelNode,
    labels,
    plotHeight,
    rowHeight,
    alignToTop = false,
    offsetTop = 0
  ) => {
    rowLabelNode.innerHTML = "";
    rowLabelNode.style.gridTemplateRows = `repeat(${labels.length}, ${rowHeight}px)`;
    rowLabelNode.style.paddingTop = alignToTop
      ? `${Math.max(0, offsetTop)}px`
      : `${Math.max(0, offsetTop + (plotHeight - labels.length * rowHeight) / 2)}px`;

    labels.forEach((label) => {
      const item = document.createElement("div");
      item.className = "sample-distance-row-label";
      item.style.height = `${rowHeight}px`;
      item.title = label;

      const text = document.createElement("span");
      text.className = "sample-distance-row-sample";
      text.textContent = label;
      item.append(text);

      rowLabelNode.append(item);
    });
  };

  const computeHeatmapSize = (availableWidth, sampleCount) => {
    const plotWidth = Math.max(
      180,
      availableWidth - heatmapPadding.left - heatmapPadding.right
    );
    const cellWidth = plotWidth / Math.max(sampleCount, 1);
    const rowHeight = clamp(Math.round(cellWidth * 0.76), 28, 52);
    return {
      width: plotWidth,
      height: Math.max(240, sampleCount * rowHeight),
      rowHeight,
    };
  };

  const computePcaSize = (availableWidth, padding) => {
    const width = Math.max(220, availableWidth - padding.left - padding.right);
    return {
      width,
      height: clamp(Math.round(width * 0.62), 240, 360),
    };
  };

  const computeCountMatrixSize = (availableWidth, sampleCount, featureCount) => {
    const plotWidth = Math.max(
      260,
      availableWidth - countMatrixPadding.left - countMatrixPadding.right
    );
    const rowHeight = clamp(Math.round(640 / Math.max(featureCount, 1)), 8, 18);
    return {
      width: plotWidth,
      height: featureCount * rowHeight,
      rowHeight,
    };
  };

  const initializeHeatmap = async () => {
    const container = document.getElementById("sample-distance-view");
    const plotShell = container ? container.parentElement : null;
    const rowLabelNode = document.getElementById("sample-distance-row-labels");
    const legendNode = document.getElementById("sample-distance-legend");
    const specNode = document.getElementById("vega_sample_distance_spec");
    const dataPathNode = document.getElementById("sample_distances_data_path");
    const orderNode = document.getElementById("sample_distance_order");
    const errorNode = document.getElementById("sample-distance-error");

    if (
      !container || !plotShell || !rowLabelNode || !legendNode || !specNode
      || !dataPathNode || !orderNode
    ) {
      return;
    }

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

      renderDistanceLegend(legendNode, cells);

      let embedResult = null;
      let currentWidth = 0;
      let currentHeight = 0;

      const applySpecConfig = (renderSpec, plotSize) => {
        renderSpec.autosize = "none";
        renderSpec.width = plotSize.width;
        renderSpec.height = plotSize.height;
        renderSpec.padding = heatmapPadding;
        renderSpec.data[0].values = cells;
        renderSpec.legends = [];
        renderSpec.axes = Array.isArray(renderSpec.axes)
          ? renderSpec.axes.filter((axis) => axis.scale === "x")
          : [];

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
          xAxis.labelFontSize = sampleOrder.length >= 20 ? 11 : 13;
          xAxis.labelAngle = -90;
          xAxis.labelAlign = "right";
          xAxis.labelBaseline = "middle";
          xAxis.labelPadding = 8;
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
        const availableWidth = measureWidth(plotShell);
        if (availableWidth < 220) {
          return;
        }

        const plotSize = computeHeatmapSize(availableWidth, sampleOrder.length);
        renderRowLabels(
          rowLabelNode,
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

      scheduleViewResize(async () => {
        try {
          await renderOrResize();
        } catch (error) {
          showError(errorNode, "The sample-distance heatmap", error);
        }
      }, [plotShell, plotShell.parentElement]);
    } catch (error) {
      showError(errorNode, "The sample-distance heatmap", error);
    }
  };

  const initializePca = async () => {
    const container = document.getElementById("sample-pca-view");
    const plotShell = container ? container.parentElement : null;
    const specNode = document.getElementById("vega_sample_pca_spec");
    const dataPathNode = document.getElementById("sample_pca_data_path");
    const errorNode = document.getElementById("sample-pca-error");

    if (!container || !plotShell || !specNode || !dataPathNode) {
      return;
    }

    try {
      const baseSpec = parseJsonNode(specNode, null);
      const dataPath = parseJsonNode(dataPathNode, "");
      if (!baseSpec || !dataPath) {
        return;
      }

      const response = await fetch(dataPath);
      if (!response.ok) {
        throw new Error(`Failed to load ${dataPath}: ${response.status}`);
      }

      const payload = await response.json();
      let points = Array.isArray(payload.points) ? payload.points : [];
      if (points.length === 0) {
        return;
      }

      const percentVariance = payload.percent_variance || {};
      const groupField = typeof payload.group_field === "string"
        ? payload.group_field
        : "";
      const groupLabel = typeof payload.group_label === "string"
        ? payload.group_label
        : "";
      const distinctGroups = Array.from(
        new Set(points.map((point) => point.group_value).filter(Boolean))
      );
      const hasGrouping = groupLabel !== "" && distinctGroups.length > 1;
      const groupPalette = hasGrouping
        ? buildGroupPalette(points, groupField)
        : { domain: [], colors: new Map() };
      if (hasGrouping) {
        points = points.map((point) => {
          const colors = groupPalette.colors.get(point.group_value) || {};
          return {
            ...point,
            group_fill_color: colors.background || "#dbeafe",
            group_stroke_color: colors.border || "#93c5fd",
          };
        });
      }

      let embedResult = null;
      let currentWidth = 0;
      let currentHeight = 0;
      const pcaPadding = {
        left: 60,
        right: 12,
        top: hasGrouping ? 60 : 12,
        bottom: 50,
      };

      const applySpecConfig = (renderSpec, plotSize) => {
        renderSpec.autosize = "none";
        renderSpec.width = plotSize.width;
        renderSpec.height = plotSize.height;
        renderSpec.data[0].values = points;
        renderSpec.padding = pcaPadding;

        const xAxis = renderSpec.axes.find((axis) => axis.scale === "x");
        const yAxis = renderSpec.axes.find((axis) => axis.scale === "y");
        if (xAxis) {
          xAxis.labelFontSize = 13;
          xAxis.titleFontSize = 15;
          xAxis.title = `PC1: ${formatPercent(Number(percentVariance.PC1))}% variance`;
          xAxis.labelPadding = 4;
          xAxis.titlePadding = 6;
        }
        if (yAxis) {
          yAxis.labelFontSize = 13;
          yAxis.titleFontSize = 15;
          yAxis.title = `PC2: ${formatPercent(Number(percentVariance.PC2))}% variance`;
          yAxis.labelPadding = 4;
          yAxis.titlePadding = 4;
        }

        if (hasGrouping) {
          const colorScale = renderSpec.scales.find((scale) => scale.name === "color");
          if (colorScale) {
            colorScale.domain = groupPalette.domain;
            colorScale.range = groupPalette.domain.map((value) => (
              groupPalette.colors.get(value)?.background || "#dbeafe"
            ));
          }
          if (Array.isArray(renderSpec.legends) && renderSpec.legends[0]) {
            renderSpec.legends[0].title = groupLabel;
            renderSpec.legends[0].labelFontSize = 13;
            renderSpec.legends[0].titleFontSize = 15;
          }
          if (renderSpec.marks?.[0]?.encode?.update?.fill) {
            renderSpec.marks[0].encode.update.fill = {
              scale: "color",
              field: "group_value",
            };
          }
          if (renderSpec.marks?.[0]?.encode?.update) {
            renderSpec.marks[0].encode.update.stroke = {
              field: "group_stroke_color",
            };
            renderSpec.marks[0].encode.update.strokeWidth = {
              value: 1.5,
            };
          }
        } else {
          renderSpec.legends = [];
          if (renderSpec.marks?.[0]?.encode?.update?.fill) {
            renderSpec.marks[0].encode.update.fill = { value: "#2563eb" };
          }
          if (renderSpec.marks?.[0]?.encode?.update) {
            renderSpec.marks[0].encode.update.stroke = { value: "#ffffff" };
            renderSpec.marks[0].encode.update.strokeWidth = { value: 1.5 };
          }
        }
      };

      const renderOrResize = async () => {
        const availableWidth = measureWidth(plotShell);
        if (availableWidth < 220) {
          return;
        }

        const plotSize = computePcaSize(availableWidth, pcaPadding);
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
          embedResult = await vegaEmbed("#sample-pca-view", renderSpec, {
            renderer: "canvas",
            actions: {
              export: { png: true, svg: true },
              source: false,
              compiled: false,
              editor: false,
            },
            downloadFileName: "deseq2-sample-pca-plot",
            scaleFactor: { png: 2, svg: 1 },
          });
          return;
        }

        embedResult.view.width(currentWidth);
        embedResult.view.height(currentHeight);
        await embedResult.view.runAsync();
      };

      scheduleViewResize(async () => {
        try {
          await renderOrResize();
        } catch (error) {
          showError(errorNode, "The sample PCA plot", error);
        }
      }, [plotShell, plotShell.parentElement]);
    } catch (error) {
      showError(errorNode, "The sample PCA plot", error);
    }
  };

  const initializeCountMatrixHeatmap = async () => {
    const container = document.getElementById("count-matrix-heatmap-view");
    const plotShell = container ? container.parentElement : null;
    const sampleLabelNode = document.getElementById("count-matrix-sample-labels");
    const legendNode = document.getElementById("count-matrix-heatmap-legend");
    const specNode = document.getElementById("vega_count_matrix_heatmap_spec");
    const dataPathNode = document.getElementById("count_matrix_heatmap_data_path");
    const errorNode = document.getElementById("count-matrix-heatmap-error");

    if (
      !container || !plotShell || !sampleLabelNode || !legendNode || !specNode
      || !dataPathNode
    ) {
      return;
    }

    try {
      const baseSpec = parseJsonNode(specNode, null);
      const dataPath = parseJsonNode(dataPathNode, "");
      if (!baseSpec || !dataPath) {
        return;
      }

      const response = await fetch(dataPath);
      if (!response.ok) {
        throw new Error(`Failed to load ${dataPath}: ${response.status}`);
      }

      const payload = await response.json();
      const cells = Array.isArray(payload.cells) ? payload.cells : [];
      const sampleOrder = Array.isArray(payload.sample_order) ? payload.sample_order : [];
      const featureOrder = Array.isArray(payload.feature_order) ? payload.feature_order : [];
      const samples = Array.isArray(payload.samples) ? payload.samples : [];
      if (cells.length === 0 || sampleOrder.length === 0 || featureOrder.length === 0) {
        return;
      }

      const parsedSampleRows = (
        samples.length === sampleOrder.length
          ? samples
          : sampleOrder.map((sampleId) => ({
              sample_id: sampleId,
              sample_label: sampleId,
              sample_metadata: "",
            }))
      ).map((sample) => ({
        metadata: parseMetadataText(sample.sample_metadata),
      }));
      const metadataColors = buildMetadataColorMap(parsedSampleRows);

      renderHeatmapLegend(
        legendNode,
        cells.map((cell) => cell.value),
        "VST count"
      );

      let embedResult = null;
      let currentWidth = 0;
      let currentHeight = 0;

      const applySpecConfig = (renderSpec, plotSize) => {
        renderSpec.width = plotSize.width;
        renderSpec.height = plotSize.height;
        renderSpec.padding = countMatrixPadding;
        renderSpec.data[0].values = cells;
        renderSpec.legends = [];
        renderSpec.axes = [];

        const xScale = renderSpec.scales.find((scale) => scale.name === "x");
        const yScale = renderSpec.scales.find((scale) => scale.name === "y");
        const colorScale = renderSpec.scales.find((scale) => scale.name === "color");
        if (xScale) {
          xScale.domain = sampleOrder;
          xScale.paddingInner = 0;
          xScale.paddingOuter = 0;
        }
        if (yScale) {
          yScale.domain = featureOrder;
          yScale.paddingInner = 0;
          yScale.paddingOuter = 0;
        }
        if (colorScale) {
          colorScale.zero = false;
        }

        if (renderSpec.marks?.[0]?.encode?.update?.x) {
          renderSpec.marks[0].encode.update.x.field = "sample_id";
        }
        if (renderSpec.marks?.[0]?.encode?.update?.y) {
          renderSpec.marks[0].encode.update.y.field = "feature_id";
        }
        if (renderSpec.marks?.[0]?.encode?.update?.fill) {
          renderSpec.marks[0].encode.update.fill = {
            scale: "color",
            field: "value",
          };
        }
        if (renderSpec.marks?.[0]?.encode?.enter) {
          renderSpec.marks[0].encode.enter.tooltip = {
            signal: "{'Feature': datum.feature_id, 'Sample': datum.sample_id, 'Sample metadata': datum.sample_metadata || 'NA', 'VST count': datum.value == null ? 'NA' : format(datum.value, '.4f')}"
          };
        }
      };

      const renderOrResize = async () => {
        const availableWidth = measureWidth(plotShell);
        if (availableWidth < 260) {
          return;
        }

        const plotSize = computeCountMatrixSize(
          availableWidth,
          sampleOrder.length,
          featureOrder.length
        );
        renderCountMatrixSampleLabels(
          sampleLabelNode,
          samples.length === sampleOrder.length
            ? samples
            : sampleOrder.map((sampleId) => ({
                sample_id: sampleId,
                sample_label: sampleId,
                sample_metadata: "",
              })),
          metadataColors,
          plotSize.width
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
          embedResult = await vegaEmbed("#count-matrix-heatmap-view", renderSpec, {
            renderer: "canvas",
            actions: {
              export: { png: true, svg: true },
              source: false,
              compiled: false,
              editor: false,
            },
            downloadFileName: "deseq2-count-matrix-heatmap",
            scaleFactor: { png: 2, svg: 1 },
          });
          return;
        }

        embedResult.view.width(currentWidth);
        embedResult.view.height(currentHeight);
        await embedResult.view.runAsync();
      };

      scheduleViewResize(async () => {
        try {
          await renderOrResize();
        } catch (error) {
          showError(errorNode, "The count-matrix heatmap", error);
        }
      }, [plotShell, plotShell.parentElement]);
    } catch (error) {
      showError(errorNode, "The count-matrix heatmap", error);
    }
  };

  await Promise.allSettled([
    initializeHeatmap(),
    initializePca(),
    initializeCountMatrixHeatmap(),
  ]);
});
