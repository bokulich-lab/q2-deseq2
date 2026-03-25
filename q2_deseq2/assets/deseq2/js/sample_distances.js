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
    const cellSize = measureCellSize(sampleOrder.length);
    const heatmapSide = Math.max(320, sampleOrder.length * cellSize);

    spec.width = heatmapSide;
    spec.height = heatmapSide;
    spec.data[0].values = cells;

    const xScale = spec.scales.find((scale) => scale.name === "x");
    const yScale = spec.scales.find((scale) => scale.name === "y");
    if (xScale) {
      xScale.domain = sampleOrder;
    }
    if (yScale) {
      yScale.domain = sampleOrder;
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
