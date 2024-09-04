import dash_bootstrap_components as dbc

NavBar = dbc.Nav(
    [
        dbc.NavLink("Home", href="/home", style={'font-weight': 'bold', 'color': 'white'}),
        dbc.NavLink("Advance Tools", href="/advance-tools", style={'font-weight': 'bold', 'color': 'white'}),
        dbc.NavLink("Research Assistant", href="/research-assistant", style={'font-weight': 'bold', 'color': 'white'}),
        dbc.NavLink("About", href="/about", style={'font-weight': 'bold', 'color': 'white'}),
    ],
    style={
        'backgroundColor': 'black',
        'width': '100%',
        'padding': '10px',
    }
)
