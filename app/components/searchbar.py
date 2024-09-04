import dash_bootstrap_components as dbc

# Define the search bar
search_bar = dbc.Row(
    [
        dbc.Col(dbc.Input(id="search-input", type="search", placeholder="Search")),  # Correct ID here
        dbc.Col(
            dbc.Button(
                "Search", color="primary", className="ms-2", n_clicks=0, id="search-button"  # Correct ID here
            ),
            width="auto",
        ),
    ],
    className="g-0 ms-auto flex-nowrap mt-3 mt-md-0",
    align="center",
)

# Define a toast component
search_toast = dbc.Toast(
    "Please enter text to search.",
    id="search-toast",  # Correct ID here
    header="Input Error",
    icon="danger",
    duration=4000,  # Duration in milliseconds (4 seconds)
    is_open=False,  # Initially closed
    style={"position": "fixed", "top": 10, "right": 10, "width": 350},
)
