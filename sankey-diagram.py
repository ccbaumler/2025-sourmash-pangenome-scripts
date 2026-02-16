#! /usr/bin/env python

import polars as pl
from polars import col
import plotly.graph_objects as go
import plotly.io as py

# Load the CSV file into a DataFrame
df = pl.read_csv('dbs/gtdb-rs214.lineages.csv')

search_term = 'Fusobacterium'

header = df.columns

filtered_rows = df.filter(col('species').str.contains(search_term))

# Create lists to hold the source, target, and value data for the Sankey diagram
sources = [] # Each index contains one source to be matched to a target in the same index
targets = [] # Contains all the targets for the sankey graph
values = [] # All 1s since each ident is unique

# Define the hierarchy
hierarchy = ['family', 'genus', 'species', 'ident']

# Create a labels array to store all the different names for each level
labels = []
# Iterate over each row in the DataFrame to populate the sources, targets, and values lists
# Loop through each row
for row in filtered_rows.iter_rows():
    # Loop each family -> genus -> species in the row
    for i, level in enumerate(hierarchy[:-1]):
        # Set the current and next level according to the hierarchy
        current_label = row[i]  # Access column value by index
        next_label = row[i+1]   # Access next column value by index

        # Add labels if does not exist
        if current_label not in labels:
            labels.append(current_label)

        if next_label not in labels:
            labels.append(next_label)

        # Set the source to the current label
        sources.append(labels.index(current_label))
        # Set the target to the next label
        targets.append(labels.index(next_label))
        # Since each ident is unique there is only 1
        values.append(1)

# Create a Sankey diagram figure using Plotly
fig = go.Figure(data=[go.Sankey(
    # Set details for all the nodes (labels)
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=labels  # Use the updated labels list
    ),
    # Use the generated source target and value array for the graph
    link=dict(
        source=sources,
        target=targets,
        value=values
    )
)])

# Adjust the height of the graph to scale with the amount of labels
height_graph = len(labels) * 10

# Update layout
fig.update_layout(title_text="Sankey Diagram", font_size=10, height=height_graph)

# Save the Sankey diagram as an HTML file
py.write_html(fig, file='sankey_diagram.html', auto_open=True)
