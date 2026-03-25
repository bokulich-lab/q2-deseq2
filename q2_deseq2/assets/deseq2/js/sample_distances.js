document.addEventListener("DOMContentLoaded", async () => {
  const container = document.getElementById("sample-distance-view");
  const specNode = document.getElementById("vega_sample_distance_spec");
  const dataPathNode = document.getElementById("sample_distances_data_path");
  const annotationPathNode = document.getElementById("sample_annotations_data_path");
  const orderNode = document.getElementById("sample_distance_order");
  const labelNode = document.getElementById("sample-distance-annotation-labels");
  const spacerNode = document.getElementById("sample-distance-annotation-spacer");
  const topAnnotationsNode = document.getElementById("sample-distance-top-annotations");
  const bottomSpacerNode = document.getElementById("sample-distance-bottom-spacer");
  const leftAnnotationsNode = document.getElementById("sample-distance-left-annotations");
  const legendNode = document.getElementById("sample-distance-annotation-legend");
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

  const annotationStripSize = 18;
  const annotationPalette = [
    "#0f766e",
    "#2563eb",
    "#dc2626",
    "#7c3aed",
    "#ea580c",
    "#059669",
    "#d97706",
    "#be123c",
    "#4f46e5",
    "#0f766e",
  ];

  const toggleAnnotationNodes = (hasAnnotations) => {
    [
      labelNode,
      spacerNode,
      topAnnotationsNode,
      bottomSpacerNode,
      leftAnnotationsNode,
      legendNode,
    ].forEach((node) => {
      if (node) {
        node.hidden = !hasAnnotations;
      }
    });
  };

  const buildAnnotationMaps = (fields, records) => {
    const annotationsBySample = new Map();
    const valuesByField = new Map(fields.map((field) => [field.field, []]));
    const seenValuesByField = new Map(fields.map((field) => [field.field, new Set()]));

    records.forEach((record) => {
      const field = typeof record.field === "string" ? record.field : "";
      const sampleId = typeof record.sample_id === "string" ? record.sample_id : "";
      if (!field || !sampleId) {
        return;
      }

      const rawValue = record.value == null ? "" : String(record.value).trim();
      const value = rawValue || "NA";

      if (!annotationsBySample.has(sampleId)) {
        annotationsBySample.set(sampleId, {});
      }
      annotationsBySample.get(sampleId)[field] = value;

      if (!valuesByField.has(field)) {
        valuesByField.set(field, []);
        seenValuesByField.set(field, new Set());
      }
      if (!seenValuesByField.get(field).has(value)) {
        valuesByField.get(field).push(value);
        seenValuesByField.get(field).add(value);
      }
    });

    return { annotationsBySample, valuesByField };
  };

  const buildColorLookup = (fields, valuesByField) => {
    const colorByKey = new Map();

    fields.forEach((field, fieldIndex) => {
      const values = valuesByField.get(field.field) || [];
      values.forEach((value, valueIndex) => {
        const color = value === "NA"
          ? "#e5e7eb"
          : annotationPalette[(fieldIndex * 3 + valueIndex) % annotationPalette.length];
        colorByKey.set(`${field.field}::${value}`, color);
      });
    });

    return colorByKey;
  };

  const formatSampleAnnotations = (sampleId, fields, annotationsBySample) => {
    const sampleAnnotations = annotationsBySample.get(sampleId) || {};
    return fields.map((field) => (
      `${field.label || field.field}: ${sampleAnnotations[field.field] || "NA"}`
    )).join("; ");
  };

  const renderAnnotationLabels = (fields) => {
    if (!labelNode) {
      return;
    }

    labelNode.innerHTML = "";
    const list = document.createElement("div");
    list.className = "sample-distance-annotation-label-list";
    list.style.gridTemplateRows = `repeat(${fields.length}, ${annotationStripSize}px)`;

    fields.forEach((field) => {
      const label = document.createElement("div");
      label.className = "sample-distance-annotation-label";
      label.style.height = `${annotationStripSize}px`;
      label.textContent = field.label || field.field;
      list.append(label);
    });

    labelNode.append(list);
  };

  const renderTopAnnotations = (
    fields,
    sampleOrder,
    cellSize,
    annotationsBySample,
    colorByKey
  ) => {
    if (!topAnnotationsNode || !spacerNode) {
      return;
    }

    topAnnotationsNode.innerHTML = "";
    spacerNode.style.width = `${fields.length * annotationStripSize}px`;

    fields.forEach((field) => {
      const row = document.createElement("div");
      row.className = "sample-distance-annotation-row";
      row.style.gridTemplateColumns = `repeat(${sampleOrder.length}, ${cellSize}px)`;

      sampleOrder.forEach((sampleId) => {
        const value = (annotationsBySample.get(sampleId) || {})[field.field] || "NA";
        const swatch = document.createElement("div");
        swatch.className = "sample-distance-annotation-swatch";
        if (value === "NA") {
          swatch.classList.add("is-missing");
        }
        swatch.style.background = colorByKey.get(`${field.field}::${value}`) || "#e5e7eb";
        swatch.style.height = `${annotationStripSize}px`;
        swatch.title = `${sampleId}\n${field.label || field.field}: ${value}`;
        row.append(swatch);
      });

      topAnnotationsNode.append(row);
    });
  };

  const renderLeftAnnotations = (
    fields,
    sampleOrder,
    cellSize,
    annotationsBySample,
    colorByKey
  ) => {
    if (!leftAnnotationsNode) {
      return;
    }

    leftAnnotationsNode.innerHTML = "";
    leftAnnotationsNode.style.gridTemplateColumns = `repeat(${fields.length}, ${annotationStripSize}px)`;
    leftAnnotationsNode.style.gridTemplateRows = `repeat(${sampleOrder.length}, ${cellSize}px)`;

    sampleOrder.forEach((sampleId) => {
      fields.forEach((field) => {
        const value = (annotationsBySample.get(sampleId) || {})[field.field] || "NA";
        const swatch = document.createElement("div");
        swatch.className = "sample-distance-annotation-swatch";
        if (value === "NA") {
          swatch.classList.add("is-missing");
        }
        swatch.style.background = colorByKey.get(`${field.field}::${value}`) || "#e5e7eb";
        swatch.style.width = `${annotationStripSize}px`;
        swatch.style.height = `${cellSize}px`;
        swatch.title = `${sampleId}\n${field.label || field.field}: ${value}`;
        leftAnnotationsNode.append(swatch);
      });
    });
  };

  const renderAnnotationLegend = (fields, valuesByField, colorByKey) => {
    if (!legendNode) {
      return;
    }

    legendNode.innerHTML = "";
    fields.forEach((field) => {
      const values = valuesByField.get(field.field) || [];
      if (values.length === 0) {
        return;
      }

      const group = document.createElement("div");
      group.className = "sample-distance-annotation-legend-group";

      const title = document.createElement("span");
      title.className = "sample-distance-annotation-legend-title";
      title.textContent = field.label || field.field;
      group.append(title);

      const items = document.createElement("div");
      items.className = "sample-distance-annotation-legend-items";

      values.forEach((value) => {
        const item = document.createElement("div");
        item.className = "sample-distance-annotation-legend-item";

        const swatch = document.createElement("span");
        swatch.className = "sample-distance-annotation-legend-swatch";
        swatch.style.background = colorByKey.get(`${field.field}::${value}`) || "#e5e7eb";
        item.append(swatch);

        const label = document.createElement("span");
        label.textContent = value;
        item.append(label);
        items.append(item);
      });

      group.append(items);
      legendNode.append(group);
    });
  };

  try {
    const spec = parseJsonNode(specNode, null);
    const sampleOrder = parseJsonNode(orderNode, []);
    const dataPath = parseJsonNode(dataPathNode, "");
    const annotationDataPath = annotationPathNode
      ? parseJsonNode(annotationPathNode, "")
      : "";
    if (!spec || !Array.isArray(sampleOrder) || sampleOrder.length === 0 || !dataPath) {
      return;
    }

    const response = await fetch(dataPath);
    if (!response.ok) {
      throw new Error(`Failed to load ${dataPath}: ${response.status}`);
    }

    const cells = await response.json();
    let annotationPayload = { fields: [], records: [] };
    if (annotationDataPath) {
      const annotationResponse = await fetch(annotationDataPath);
      if (annotationResponse.ok) {
        annotationPayload = await annotationResponse.json();
      }
    }

    const fields = Array.isArray(annotationPayload.fields)
      ? annotationPayload.fields
        .filter((field) => field && typeof field.field === "string")
        .map((field) => ({
          ...field,
          label: typeof field.label === "string" && field.label.trim()
            ? field.label
            : field.field,
        }))
      : [];
    const records = Array.isArray(annotationPayload.records)
      ? annotationPayload.records
      : [];
    const hasAnnotations = fields.length > 0 && records.length > 0;
    toggleAnnotationNodes(hasAnnotations);

    const { annotationsBySample, valuesByField } = buildAnnotationMaps(fields, records);
    const colorByKey = buildColorLookup(fields, valuesByField);
    const cellSize = measureCellSize(sampleOrder.length);
    const heatmapSide = Math.max(320, sampleOrder.length * cellSize);
    container.style.width = `${heatmapSide}px`;
    container.style.minHeight = `${heatmapSide}px`;

    if (hasAnnotations) {
      renderAnnotationLabels(fields);
      renderTopAnnotations(
        fields,
        sampleOrder,
        cellSize,
        annotationsBySample,
        colorByKey
      );
      renderLeftAnnotations(
        fields,
        sampleOrder,
        cellSize,
        annotationsBySample,
        colorByKey
      );
      renderAnnotationLegend(fields, valuesByField, colorByKey);
    }

    const enrichedCells = cells.map((cell) => ({
      ...cell,
      sample_x_metadata: hasAnnotations
        ? formatSampleAnnotations(cell.sample_x, fields, annotationsBySample)
        : "",
      sample_y_metadata: hasAnnotations
        ? formatSampleAnnotations(cell.sample_y, fields, annotationsBySample)
        : "",
    }));

    spec.width = heatmapSide;
    spec.height = heatmapSide;
    spec.data[0].values = enrichedCells;

    const xScale = spec.scales.find((scale) => scale.name === "x");
    const yScale = spec.scales.find((scale) => scale.name === "y");
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

    if (spec.marks?.[0]?.encode?.enter) {
      spec.marks[0].encode.enter.tooltip = hasAnnotations
        ? {
            signal: "{'Sample A': datum.sample_x, 'Sample A metadata': datum.sample_x_metadata, 'Sample B': datum.sample_y, 'Sample B metadata': datum.sample_y_metadata, 'Distance': datum.distance == null ? 'NA' : format(datum.distance, '.4f')}"
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
