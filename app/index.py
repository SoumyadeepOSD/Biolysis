from pages.research_assistant import ResearchAssistantPage
from dash import dcc, html, Input, Output, State
from pages.advance_tools import AdvanceToolsPage
from utils.extract import generate_output
import dash_bootstrap_components as dbc
from components.navbar import NavBar
from pages.about import AboutPage
from pages.home import HomePage
import plotly.graph_objs as go
import dash

class MainApplication:
    def __init__(self):
        self.__app = dash.Dash(
            __name__,
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            suppress_callback_exceptions=True,
            use_pages=True
        )
        self.set_layout()

    @property
    def app(self):
        return self.__app

    def set_layout(self):
        self.app.layout = dbc.Container(
            children=[
                html.Div(style={'height': '10px'}),
                NavBar,  # Ensure NavBar is callable
                html.Div(style={'height': '10px'}),
                dcc.Location(id='url', refresh=False),  # Location component for URL routing
                html.Div(id='page-content')  # Placeholder for rendering page content
            ],
            style={
                'backgroundColor': 'white',
                'width': '100vw',
                'height': '100vh',
            }
        )

        # Callback to handle page routing
        @self.app.callback(
            Output('page-content', 'children'),
            [Input('url', 'pathname')]
        )
        def display_page(pathname):
            if pathname == '/':
                return HomePage()
            elif pathname == '/advance-tools':
                return AdvanceToolsPage()
            elif pathname == '/research-assistant':
                return ResearchAssistantPage()
            elif pathname == '/about':
                return AboutPage()
            else:
                return html.Div("404 Page Not Found", style={'textAlign': 'center'})

        # Callback to handle toast visibility and update search results
        @self.app.callback(
            [
                Output('search-toast', 'is_open'),
                Output('search-results', 'children'),
                Output('2d-drawing', 'children'),
                Output('3d-drawing', 'children')
            ],
            [Input('search-button', 'n_clicks')],
            [State('search-input', 'value')]
        )
        def SearchResults(n_clicks, search_value):
            if n_clicks > 0:
                if not search_value:
                    return True, None, None, None  # Show toast and no search results
                else:
                    is_open = False
                    # Get the search results
                    compound_info, img_element, graph = generate_output(n_clicks, search_value)
                    return is_open, compound_info, img_element, graph
            return False, dbc.Alert("Search by compound name only e.g., Benzene", color="info"), None, None


# Initialize the application
Application = MainApplication()
app = Application.app.server

if __name__ == "__main__":
    Application.app.run_server(debug=True)
