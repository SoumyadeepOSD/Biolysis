from phi.assistant import Assistant
from phi.llm.groq import Groq
import os
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from phi.tools.wikipedia import WikipediaTools
from phi.tools.tavily import TavilyTools
from phi.tools.website import WebsiteTools
from pydantic import BaseModel, Field
from dotenv import load_dotenv
from phi.llm.google import Gemini
import re   # Import re module for regular expressions

load_dotenv()

os.environ["GROQ_API_KEY"] = os.getenv("GROQ_API_KEY")
os.environ["TAVILY_API_KEY"] = os.getenv("TAVILY_API_KEY")
response = ""

st.header("Autonomous Assistant for MicroBiologist")

compound = ["C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O"]
ms = [Chem.MolFromSmiles(ms) for ms in compound]
img = Draw.MolsToGridImage(ms, molsPerRow=1)
st.write(img)

class CompoundStructure(BaseModel):
    name: str = Field(..., description="Proper market name of the compound")
    iupacName: str = Field(..., description="IUPAC name of the compound")
    genre: str = Field(..., description="Characteristics of compound")
    molecular_structure: str = Field(..., description="Provide the molecular structure in canonical SMILES format only, which can be used directly with RDKit")
    description: str = Field(..., description="3-sentence description for the compound. Make it informative")

models = [
    "llama3-8b-8192",
    "mixtral-8x7b-32768",
    "gemini-1.5-flash"
]

microbiologist = Assistant(
    llm=Groq(model=models[1]),
    description="You help scientists to identify and diagnose molecular structures in canonical SMILES format.",
    tools=[
            WikipediaTools(), 
            TavilyTools(),
        ],
    instructions=["You are providing details of organic compounds in JSON format only with the following fields given in CompoundStructure. Ensure the molecular_structure is in canonical SMILES format."],
    debug_mode=True,
    show_tool_calls=True,
    output_model=CompoundStructure
)

secratary = Assistant(
    llm=Groq(model=models[1]),
    description="You extract the canonical smile of any compound, name given by user",
    tools=[WebsiteTools(),TavilyTools()],
    instructions=["You go to https://pubchem.ncbi.nlm.nih.gov/compound and tell me the canonical smile name of compound and print into console"],
    debug_mode=True,
    show_tool_calls=True,
)

secratary.print_response("Give me canonical form of Doxorubicin", markdown=True)

# prompt = st.text_input("Please Enter your query")
# if prompt != "":
#     button = st.button("Submit")

# if prompt is not None and prompt != "" and button:
#     detailed_prompt = f"Provide the details of the compound '{prompt}' including its molecular structure in canonical SMILES format. Ensure that the JSON includes the following fields: name, iupacName, genre, molecular_structure (in canonical SMILES format), and description."
#     response = ""
#     for delta in microbiologist.run(detailed_prompt):
#         response += delta

#     # Extract the JSON portion of the response
#     start_index = response.find("{")
#     end_index = response.rfind("}") + 1

#     if response != "":
#         if start_index != -1 and end_index != -1:
#             json_response = response[start_index:end_index]
            
#             # Using regular expression to extract key-value pairs
#             pattern = r'\"(.*?)\":\s*\"(.*?)\"'
#             matches = re.findall(pattern, json_response)
            
#             # Convert matches to a dictionary
#             extracted_data = {key: value for key, value in matches}
            
#             # Display the extracted data in Streamlit
#             st.write("Extracted Key-Value Pairs:")
#             for key, value in extracted_data.items():
#                 st.write(f"{key}: {value}")

#             # Extract SMILES structure
#             smiles_structure = extracted_data.get("molecular_structure")
#             if smiles_structure:
#                 # Convert SMILES to molecule and draw
#                 mol = Chem.MolFromSmiles(smiles_structure)
#                 if mol:
#                     img = Draw.MolToImage(mol)
#                     st.image(img, caption="Molecular Structure")
#                 else:
#                     st.error("The provided molecular structure is not a valid SMILES.")
#             else:
#                 st.error("No molecular structure found in the response.")
#         else:
#             st.error("No valid JSON data found in the response.")
