import dash_bootstrap_components as dbc
from dash import html, dcc
from components.searchbar import search_bar, search_toast

def HomePage():
    return dbc.Card(
        children=[
            search_bar,
            search_toast,
            html.Div(id='compound-info', style={'marginTop': '20px'}),  # Placeholder for compound info
            dbc.Row(
                id='drawings-container',
                children=[
                    dbc.Col(html.Div(id='2d-drawing', style={'textAlign': 'center'}), width=12, md=6),
                    dbc.Col(html.Div(id='3d-drawing', style={'textAlign': 'center'}), width=12, md=6)
                ],
                className="g-3",  # Add some gap between the columns
                justify="center",
                align="center",
                style={'marginTop': '20px'}
            ),
            html.Div(id='search-results', style={'marginTop': '20px'})  # Placeholder for search results
        ],
        body=True,
        style={
            'border': '0.5px solid #a6a6a6', 
            'padding': '20px', 
            'height': '100%', 
            'boxShadow': '0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19)'  # Corrected shadow property
        }
    )
