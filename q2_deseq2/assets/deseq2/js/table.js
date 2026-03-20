document.addEventListener("DOMContentLoaded", async () => {
  const dataPathNode = document.getElementById("results_data_path");
  const comparisonOptionsNode = document.getElementById("comparison_options");
  const defaultComparisonNode = document.getElementById("default_comparison");
  const comparisonSelect = document.getElementById("results-comparison-select");
  const comparisonLabelNode = document.getElementById("results-selected-comparison");
  const alphaInput = document.getElementById("results-alpha-cutoff-input");
  const lfcInput = document.getElementById("results-lfc-cutoff-input");
  const resetButton = document.getElementById("results-filter-reset");
  const errorNode = document.getElementById("results-table-error");
  const tableNode = document.getElementById("results-table");

  if (
    !dataPathNode || !comparisonOptionsNode || !defaultComparisonNode
    || !comparisonSelect || !comparisonLabelNode
    || !alphaInput || !lfcInput || !tableNode
  ) {
    return;
  }

  const comparisonPersistenceKey = "q2_deseq2_selected_comparison";
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

  const readPersistedComparison = (availableComparisons) => {
    try {
      const persisted = window.localStorage.getItem(comparisonPersistenceKey);
      if (persisted && availableComparisons.includes(persisted)) {
        return persisted;
      }
    } catch {}
    return null;
  };

  const persistComparison = (comparison) => {
    try {
      window.localStorage.setItem(comparisonPersistenceKey, comparison);
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
    const configuredDefaultComparison = parseJsonNode(defaultComparisonNode, "");
    const comparisonOptions = parseJsonNode(comparisonOptionsNode, []);

    const response = await fetch(dataPath);
    if (!response.ok) {
      throw new Error(`Failed to load ${dataPath}: ${response.status}`);
    }

    const payload = await response.json();
    const columns = payload.columns || [];
    const rows = payload.data || [];
    const comparisonColumnIndex = columns.indexOf("comparison");
    const testLevelColumnIndex = columns.indexOf("test_level");
    const referenceLevelColumnIndex = columns.indexOf("reference_level");
    const padjColumnIndex = columns.indexOf("padj");
    const lfcColumnIndex = columns.indexOf("log2FoldChange");
    const rowsByComparison = new Map();

    comparisonOptions.forEach((option) => {
      rowsByComparison.set(option.comparison, []);
    });
    rows.forEach((row) => {
      const comparison = comparisonColumnIndex >= 0
        ? row[comparisonColumnIndex]
        : formatComparison(
          testLevelColumnIndex >= 0 ? row[testLevelColumnIndex] : "",
          referenceLevelColumnIndex >= 0 ? row[referenceLevelColumnIndex] : ""
        ) || configuredDefaultComparison;
      if (!rowsByComparison.has(comparison)) {
        rowsByComparison.set(comparison, []);
      }
      rowsByComparison.get(comparison).push(row);
    });

    const availableComparisons = Array.from(rowsByComparison.keys()).filter(Boolean);
    if (availableComparisons.length === 0) {
      throw new Error("No comparisons were available in the DESeq2 results table.");
    }

    availableComparisons.forEach((comparison) => {
      const option = document.createElement("option");
      option.value = comparison;
      option.textContent = comparison;
      comparisonSelect.append(option);
    });

    const hiddenColumns = new Set(["comparison", "test_level", "reference_level"]);
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

    let currentComparison = readPersistedComparison(availableComparisons)
      || (availableComparisons.includes(configuredDefaultComparison)
        ? configuredDefaultComparison
        : availableComparisons[0]);

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
      persistComparison(currentComparison);
      persistNumeric(alphaPersistenceKey, currentAlpha);
      persistNumeric(lfcPersistenceKey, currentLfc);

      const comparisonRows = rowsByComparison.get(currentComparison) || [];
      const filteredRows = comparisonRows.filter(rowPassesFilters);
      renderRows(filteredRows);

      comparisonLabelNode.textContent = currentComparison;
      comparisonSelect.value = currentComparison;
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

    comparisonSelect.addEventListener("change", () => {
      currentComparison = comparisonSelect.value;
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
