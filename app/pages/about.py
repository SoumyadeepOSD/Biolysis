import dash_bootstrap_components as dbc
from dash import html

def AboutPage():
    return dbc.Card(
        html.Div("About Us Page"),
        body=True,
        style={'border': '2px solid #007bff', 'padding': '20px', 'height': '85vh'}
    )
