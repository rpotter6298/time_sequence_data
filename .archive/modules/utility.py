import re
import pandas as pd


def text_to_table(text: str):
    rows = []
    lines = text.split("\n")
    treatment = re.search(r"Comparing treatment (.*?):", lines[0]).group(1).strip()
    for line in lines:
        metrics = re.findall(
            r"(\w+.*?) - p-value: (\d+\.\d+), power difference: ([-+]?\d*\.\d+), means: ([-+]?\d*\.\d+) vs ([-+]?\d*\.\d+)",
            line,
        )
        for metric, p_value, power_diff, mean1, mean2 in metrics:
            row = {
                "Treatment": treatment,
                "Metric": metric,
                "p-value": float(p_value),
                "Power difference": float(power_diff),
                "Mean 1": float(mean1),
                "Mean 2": float(mean2),
            }
            rows.append(row)

    df = pd.DataFrame(rows)
    return df


def text_to_table(text):
    rows = text.strip().split("\n")[1:]
    data = []
    for row in rows:
        tokens = row.split()
        treatment = tokens[0] + " " + tokens[1]
        metric = tokens[2][1:-1]
        p_value = float(tokens[-1])
        data.append([treatment, metric, p_value])
    df = pd.DataFrame(data, columns=["Treatment", "Metric", "p-value"])
    return df


text = working_model.max_growth_rate_info
text_to_table(working_model.max_growth_rate_info)
