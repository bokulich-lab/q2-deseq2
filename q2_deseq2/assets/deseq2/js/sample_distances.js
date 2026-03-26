document.addEventListener("DOMContentLoaded", async () => {
  const container = document.getElementById("sample-distance-view");
  const specNode = document.getElementById("vega_sample_distance_spec");
  const dataPathNode = document.getElementById("sample_distances_data_path");
  const orderNode = document.getElementById("sample_distance_order");
  const errorNode = document.getElementById("sample-distance-error");

  if (!container || !specNode || !dataPathNode || !orderNode) {
    return;
  }

  const parseJsonNode = (node, fallback) => {
    try {
      return JSON.parse(node.textContent);
    } catch {
      return fallback;
    }
  };

  const measureCellSize = (sampleCount) => {
    if (sampleCount >= 60) {
      return 12;
    }
    if (sampleCount >= 36) {
      return 16;
    }
    if (sampleCount >= 20) {
      return 20;
    }
    return 28;
  };

  const measureLeftPadding = (labels) => {
    const longest = labels.reduce(
      (maximum, label) => Math.max(maximum, String(label || "").length),
      0
    );
    return Math.min(620, Math.max(160, longest * 7 + 32));
  };

  try {
    const spec = parseJsonNode(specNode, null);
    const sampleOrder = parseJsonNode(orderNode, []);
    const dataPath = parseJsonNode(dataPathNode, "");
    if (!spec || !Array.isArray(sampleOrder) || sampleOrder.length === 0 || !dataPath) {
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
    const hasMetadataLabels = rowLabels.some(
      (label, index) => label !== sampleOrder[index]
    );
    const hasTooltipMetadata = cells.some((cell) => (
      (typeof cell.sample_x_metadata === "string" && cell.sample_x_metadata !== "")
      || (typeof cell.sample_y_metadata === "string" && cell.sample_y_metadata !== "")
    ));

    const cellSize = measureCellSize(sampleOrder.length);
    const heatmapSide = Math.max(320, sampleOrder.length * cellSize);
    const leftPadding = measureLeftPadding(rowLabels);

    spec.width = heatmapSide;
    spec.height = heatmapSide;
    spec.padding = {
      left: leftPadding,
      right: 60,
      top: 10,
      bottom: 130,
    };
    spec.data[0].values = cells;

    const xScale = spec.scales.find((scale) => scale.name === "x");
    const yScale = spec.scales.find((scale) => scale.name === "y");
    if (xScale) {
      xScale.domain = sampleOrder;
      xScale.paddingInner = 0;
      xScale.paddingOuter = 0;
    }
    if (yScale) {
      yScale.domain = rowLabels;
      yScale.paddingInner = 0;
      yScale.paddingOuter = 0;
    }

    const xAxis = spec.axes.find((axis) => axis.scale === "x");
    if (xAxis) {
      xAxis.labelFontSize = sampleOrder.length >= 20 ? 10 : 12;
    }

    const yAxis = spec.axes.find((axis) => axis.scale === "y");
    if (yAxis) {
      yAxis.labelFontSize = sampleOrder.length >= 20 ? 10 : 12;
      yAxis.labelLimit = Math.max(140, leftPadding - 24);
      yAxis.title = hasMetadataLabels ? "Sample / metadata" : "Sample";
    }

    if (spec.marks?.[0]?.encode?.update?.y) {
      spec.marks[0].encode.update.y.field = "sample_y_label";
    }

    if (spec.marks?.[0]?.encode?.enter) {
      spec.marks[0].encode.enter.tooltip = hasTooltipMetadata
        ? {
            signal: "{'Sample A': datum.sample_x, 'Sample A metadata': datum.sample_x_metadata || 'NA', 'Sample B': datum.sample_y, 'Sample B metadata': datum.sample_y_metadata || 'NA', 'Distance': datum.distance == null ? 'NA' : format(datum.distance, '.4f')}"
          }
        : {
            signal: "{'Sample A': datum.sample_x, 'Sample B': datum.sample_y, 'Distance': datum.distance == null ? 'NA' : format(datum.distance, '.4f')}"
          };
    }

    await vegaEmbed("#sample-distance-view", spec, {
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
  } catch (error) {
    if (errorNode) {
      errorNode.hidden = false;
      errorNode.textContent =
        "The sample-distance heatmap could not be loaded.\n\n" + String(error);
    }
  }
});
