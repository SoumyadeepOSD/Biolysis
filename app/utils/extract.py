import plotly.graph_objects as go
from rdkit.Chem import Draw
from dash import dcc, html
from rdkit import Chem
from io import BytesIO
import requests
import base64

# Atomic number to element name mapping
atomic_number_to_element = {
    1: 'Hydrogen',
    2: 'Helium',
    3: 'Lithium',
    4: 'Beryllium',
    5: 'Boron',
    6: 'Carbon',
    7: 'Nitrogen',
    8: 'Oxygen',
    9: 'Fluorine',
    10: 'Neon',
    11: 'Sodium',
    12: 'Magnesium',
    13: 'Aluminum',
    14: 'Silicon',
    15: 'Phosphorus',
    16: 'Sulfur',
    17: 'Chlorine',
    18: 'Argon',
    19: 'Potassium',
    20: 'Calcium',
    # Add other elements as needed...
}

def fetch_data(compound_name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/JSON"
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        return None

def generate_2d_image(smiles):
    try:
        compound = [smiles]
        ms = [Chem.MolFromSmiles(mol) for mol in compound]
        img = Draw.MolsToGridImage(ms, molsPerRow=1)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    except Exception as e:
        print(f"Error generating 2D image: {e}")
        return None

def generate_output(n_clicks, compound_name):
    if n_clicks > 0 and compound_name:
        json_data = fetch_data(compound_name.lower())
        if json_data:
            # Extract atomic data
            compound_info = json_data['PC_Compounds'][0]
            props = compound_info['props']

            # Extract Compound Name and SMILES
            compound_name = next((prop['value']['sval'] for prop in props if prop['urn']
                                 ['label'] == 'IUPAC Name' and prop['urn']['name'] == 'Preferred'), 'N/A')
            smiles_formula = next((prop['value']['sval'] for prop in props if prop['urn']
                                  ['label'] == 'SMILES' and prop['urn']['name'] == 'Canonical'), 'N/A')

            # Generate 2D Image
            img_str = generate_2d_image(smiles_formula)
            img_element = html.Img(src=img_str, style={'width': '300px', 'height': '300px', 'object-fit': 'contain','border': '2px solid blue', 'padding': '50px'}) if img_str else None

            # Prepare coordinates for 3D plot
            try:
                coords = compound_info['coords'][0]['conformers'][0]
                x_coords = coords['x']
                y_coords = coords['y']
                z_coords = [0.0] * len(x_coords)
                elements = compound_info['atoms']['element']
                bonds = compound_info['bonds']

                # Map atomic numbers to colors
                element_colors = {1: 'white', 6: 'black', 8: 'red'}  # H, C, O

                # Prepare data for the 3D scatter plot
                scatter_data = go.Scatter3d(
                    x=x_coords,
                    y=y_coords,
                    z=z_coords,
                    mode='markers',
                    marker=dict(
                        size=15,
                        color=[element_colors.get(ele, 'gray') for ele in elements],
                        opacity=0.8
                    ),
                    text=[f'{atomic_number_to_element.get(ele, "Unknown")} (Atomic Number: {ele})' for ele in elements],
                    customdata=elements,
                    hovertemplate='<b>%{text}</b><extra></extra>',  # Custom hover template
                )

                # Prepare data for the bonds
                bond_lines = []
                for a1, a2, order in zip(bonds['aid1'], bonds['aid2'], bonds['order']):
                    bond_lines.append(
                        go.Scatter3d(
                            x=[x_coords[a1-1], x_coords[a2-1]],
                            y=[y_coords[a1-1], y_coords[a2-1]],
                            z=[z_coords[a1-1], z_coords[a2-1]],
                            mode='lines',
                            line=dict(color='#a6a6a6', width=order*10),
                            hoverinfo='skip'
                        )
                    )

                # Combine atom and bond plots
                fig = go.Figure(
                    data=[scatter_data] + bond_lines,
                )
                fig.update_layout(
                    title=f"Compound Name: {compound_name}",
                    scene=dict(
                        xaxis_title='X',
                        yaxis_title='Y',
                        zaxis_title='Z'
                    ),
                    margin=dict(l=0, r=0, b=0, t=40)  # Remove extra margins for better border fit
                )

                # Return elements for the 2D and 3D drawings and compound info
                return (
                    html.Div([
                        html.H5(f"SMILES Formula: {smiles_formula}", style={'marginBottom': '50px'})
                    ]),
                    img_element,
                    html.Div(children=[
                        dcc.Graph(
                            figure=fig,
                            animate=True,
                            animation_options={"frame": {"redraw": True}, "transition": {"duration": 1000}},
                            style={'width': '100%', 'height': '60vh', 'border': '2px solid blue', 'marginTop': '20px'}  # Added blue border
                        )
                    ],
                    style={'width': '100%', 'height': '50%', 'display': 'flex', 'flex-direction': 'row', 'justify-content': 'center', 'align-items': 'center'}
                    )
                )
            except KeyError as e:
                print(f"Error preparing 3D plot data: {e}")
                return (
                    html.Div("Error preparing 3D plot data. Please try again."),
                    img_element,
                    None
                )
        else:
            return (
                html.Div("Error fetching data from PubChem. Please try again."),
                None,
                None
            )

    return (
        html.Div("Enter a compound name and click 'Search'."),
        None,
        None
    )
