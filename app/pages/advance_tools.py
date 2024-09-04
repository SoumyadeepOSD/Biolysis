import dash_bootstrap_components as dbc
from dash import html

def AdvanceToolsPage():
    return dbc.Card(
        html.Div("This is the Advance Tools Page"),
        body=True,
        style={'border': '2px solid #007bff', 'padding': '20px', 'height': '85vh'}
    )
