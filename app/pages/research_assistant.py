import dash_bootstrap_components as dbc
from dash import html

def ResearchAssistantPage():
    return dbc.Card(
        html.Div("This is the Research Assistant Page"),
        body=True,
        style={'border': '2px solid #007bff', 'padding': '20px', 'height': '85vh'}
    )
