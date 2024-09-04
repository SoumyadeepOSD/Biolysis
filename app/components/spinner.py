import dash_bootstrap_components as dbc
from dash import html

Spinner = html.Div(
    dbc.Spinner(color="info"),
    style={'display': 'flex', 'justify-content': 'center', 'align-items': 'center', 'height': '50%'},
)
