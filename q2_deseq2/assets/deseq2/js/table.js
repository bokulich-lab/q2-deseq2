document.addEventListener("DOMContentLoaded", async () => {
  const dataPathNode = document.getElementById("table_data_path");
  const errorNode = document.getElementById("results-table-error");
  const tableNode = document.getElementById("results-table");

  if (!dataPathNode || !tableNode) {
    return;
  }

  try {
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

    const isNumericColumn = (rows, columnIndex) => {
      return rows.every((row) => {
        const value = row[columnIndex];
        return value === null || value === "" || parseNumeric(value) !== null;
      });
    };

    const dataPath = dataPathNode.textContent.trim();
    const response = await fetch(dataPath);
    if (!response.ok) {
      throw new Error(`Failed to load ${dataPath}: ${response.status}`);
    }

    const payload = await response.json();
    const columns = payload.columns.map((name, columnIndex) => {
      const numericColumn = isNumericColumn(payload.data, columnIndex);
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

    new DataTable(tableNode, {
      columns,
      data: payload.data,
      pageLength: 25,
      scrollX: true,
    });
  } catch (error) {
    if (errorNode) {
      errorNode.hidden = false;
      errorNode.textContent =
        "The result table could not be loaded.\n\n" + String(error);
    }
  }
});
