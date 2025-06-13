from flask import Flask, render_template, jsonify, request, send_from_directory
import pandas as pd
import os

from constants import ALL_PLDDT_WINDOWS, MIN_CONTACTS_THRESHOLDS, SUPERFOLDER_TO_FASTA_AND_FOLDER

app = Flask(__name__, template_folder='web_visualization/')

@app.route('/static/merged_results/<path:filename>')
def serve_merged_results(filename):
    return send_from_directory(os.path.join('output_merged_results'), filename)

@app.route('/')
def index():
    return render_template('index.html')

# TODO: Not used for now, but might be useful in the future, will just leave here for now

# binding_domain_svg_files = []
# GAP_dfs = []
# GEF_dfs = []
# for window in ALL_PLDDT_WINDOWS:
#     if window < 0:
#         GAP_dfs.append(pd.read_csv(f"output_merged_results/GAP_merged.csv"))
#         GEF_dfs.append(pd.read_csv(f"output_merged_results/GEF_merged.csv"))
#     else:
#         GAP_dfs.append(pd.read_csv(f"output_merged_results/GAP_merged_plddt_window_{window}.csv"))
#         GEF_dfs.append(pd.read_csv(f"output_merged_results/GEF_merged_plddt_window_{window}.csv"))

# @app.route('/binding_domain_svgs')
# def binding_domain_svgs():
#     return jsonify({"binding_domain_svgs": binding_domain_svg_files})

# @app.route('/static/binding_domain_svgs/<path:filename>')
# def serve_svg(filename):
#     if filename in binding_domain_svg_files:
#         return send_from_directory(os.path.join('output_analyze_results'), filename)
    
#     return jsonify({"error": "File not found"}), 404

# @app.route('/process', methods=['POST'])
# def process_data():
#     input_val = request.json.get("input")
#     result = some_python_function(input_val)
#     return jsonify({"result": result})

# def some_python_function(val):
#     return f"Processed: {val}"

# def get_and_generate_binding_domain_svg_files():
#     # TODO: call analyze_results.py functions to generate the SVGs if they don't exist
#     for window in ALL_PLDDT_WINDOWS:
#         for contact_threshold in MIN_CONTACTS_THRESHOLDS:
#             # check if the file exists
#             for superfolder in SUPERFOLDER_TO_FASTA_AND_FOLDER.keys():
#                 file = f"{superfolder}_binding_domain_all_combinations_plddt_window_{window}_MIN_{contact_threshold}.svg"
#                 binding_domain_svg_files.append(file)
#                 if os.path.exists(os.path.join('output_analyze_results', file)):
#                     binding_domain_svg_files.append(file)
#                 else:
#                     # TODO: implement the call to generate the SVG
#                     pass

if __name__ == '__main__':
    # get_and_generate_binding_domain_svg_files()
    app.run(debug=True)
