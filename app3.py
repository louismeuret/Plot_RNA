import eventlet

eventlet.monkey_patch()
import time
from flask import (
    Flask,
    request,
    render_template,
    send_from_directory,
    redirect,
    url_for,
    session,
    send_file,
)
from werkzeug.utils import secure_filename
import uuid
import os
import glob
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.coordinates import PDB
from MDAnalysis.analysis import rms
from MDAnalysis.topology.guessers import guess_types
from MDAnalysis.analysis.align import alignto

import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
import json
import barnaba as bb
import plotly
from plotly.tools import mpl_to_plotly
from plotly.subplots import make_subplots
from flask_socketio import SocketIO, join_room, leave_room, emit
from io import BytesIO
from PIL import Image, ImageDraw
import sys
import numpy as np
import mdtraj as md
import energy_3dplot
from FoldingAnalysis.analysis import *
import plotly.graph_objects as go
import plotly.express as px
from string import ascii_uppercase
import barnaba as bb
import pickle
import io
import zipfile
from collections import defaultdict

import generate_contact_maps


app = Flask(__name__)
app.secret_key = "pi"
socketio = SocketIO(app, logger=True, engineio_logger=True)


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        directory = request.form["directory"]
        # Perform actions based on the form input, e.g., list files
        files = list_files(directory)
        return render_template("results.html", files=files)
    return render_template("index.html")


def list_files(directory):
    # Your implementation of listing files
    return os.listdir(directory)  # Simplified example




def plotly_to_json(fig):
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)


def plot_ermsd(ermsd):
    fig = go.Figure(data=go.Scatter(y=ermsd, mode="markers"))
    fig.update_layout(
        title="ERMSD Scatter Plot", xaxis_title="Frame", yaxis_title="ERMSD Value"
    )
    return fig


def plot_rmsd(rmsd):
    fig = go.Figure(data=go.Scatter(y=rmsd, mode="markers"))
    fig.update_layout(
        title="RMSD Scatter Plot", xaxis_title="Frame", yaxis_title="RMSD Value"
    )
    return fig


def plot_dotbracket(dotbracket_data):
    reverse_mapping = {}
    for frame, dotbracket in enumerate(dotbracket_data, start=1):
        reverse_mapping.setdefault(dotbracket, []).append(frame)

    structures_data = []
    for dotbracket, frames in reverse_mapping.items():
        for frame in frames:
            structures_data.append({"Frame": frame, "DotBracket": dotbracket})

    structures_df = pd.DataFrame(structures_data)
    structures_df.sort_values(by=["DotBracket", "Frame"], inplace=True)

    unique_structures = structures_df["DotBracket"].unique()
    structure_colors = [
        f"rgb{tuple(int(255 * x) for x in plt.cm.tab20b(i)[:3])}"
        for i in np.linspace(0, 1, len(unique_structures))
    ]
    color_dict = dict(zip(unique_structures, structure_colors))

    traces = []
    for i, dotbracket in enumerate(unique_structures):
        structure_df = structures_df[structures_df["DotBracket"] == dotbracket]
        x_values = []
        y_values = []
        prev_frame = None
        for _, row in structure_df.iterrows():
            frame = row["Frame"]
            if prev_frame is not None and frame != prev_frame + 1:
                x_values.append(None)
                y_values.append(None)
            x_values.append(frame)
            y_values.append(i + 1 - 0.2)
            prev_frame = frame
        trace = go.Scatter(
            x=x_values,
            y=y_values,
            mode="lines",
            line=dict(color=color_dict[dotbracket], width=8),
            name=dotbracket,
        )
        traces.append(trace)

    layout = go.Layout(
        title="Timeline of RNA Structures",
        xaxis=dict(title="Frame", showgrid=False, zeroline=False, showline=False),
        yaxis=dict(
            title="Dot-Bracket Structure",
            tickmode="array",
            tickvals=list(range(1, len(unique_structures) + 1)),
            ticktext=unique_structures,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticklen=1,
        ),
        showlegend=False,
        plot_bgcolor="rgba(255, 255, 255, 0)",
    )

    fig = go.Figure(data=traces, layout=layout)
    return fig


def plot_torsion(angles, res, torsionResidue):

    if torsionResidue.isdigit():
        residue_index = int(torsionResidue)
    else:
        residue_index = res.index(torsionResidue)

    app.logger.info(f"RESIDUE INDEX = {residue_index}")

    residue_index = 2
    specific_residue_angles = angles[:, residue_index, :]  # All frames for residue 2
    angles_names = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]

    # Define colors for each angle in hexadecimal
    colors = [
        "#636efa",
        "#ef553b",
        "#00cc96",
        "#a45bf9",
        "#fa9f5e",
        "#1ad0f2",
        "#ff6692",
    ]

    # Create subplots with 7 rows and 2 columns
    names = [
        f"Angle: {angles_names[i]} - {plot_type}"
        for i in range(7)
        for plot_type in ["Line Plot", "Histogram"]
    ]
    fig = make_subplots(
        rows=7,
        cols=2,
        shared_xaxes=False,
        subplot_titles=names,
        horizontal_spacing=0.1,  # Adjust spacing between plots
        vertical_spacing=0.1,
    )

    # Add traces for each angle (line plots in the first column and histograms in the second column)
    for i in range(7):
        y_values = specific_residue_angles[:, i]
        x_values = np.arange(len(y_values))
        valid_indices = ~np.isnan(y_values)
        color = colors[i]

        # Line plot
        fig.add_trace(
            go.Scatter(
                x=x_values[valid_indices],
                y=y_values[valid_indices],
                mode="lines",
                name=f"{angles_names[i]} - Line Plot",
                line=dict(color=color),
            ),
            row=i + 1,
            col=1,
        )

        # Histogram
        fig.add_trace(
            go.Histogram(
                x=y_values[valid_indices],
                name=f"{angles_names[i]} - Histogram",
                marker=dict(color=color),
            ),
            row=i + 1,
            col=2,
        )

    # Update layout
    fig.update_layout(
        height=1400,
        width=1600,
        title_text=f"Angles for Residue {res[residue_index]} Across All Frames",
        showlegend=False,
    )

    return fig


def plot_sec_structure(sec_structure):
    fig = go.Figure(data=go.Scatter(y=sec_structure, mode="markers"))
    fig.update_layout(
        title="Secondary Structure Plot",
        xaxis_title="Frame",
        yaxis_title="Secondary Structure",
    )
    return fig


def plot_landscapes_3D(
    energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions
):
    fig = energy_3dplot.energy_plot_3d(
        energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions
    )
    return fig

def plot_diagram_frequency(sequence, dot_bracket_list, dotbracket_native):

    def parse_dot_bracket(dot_bracket):
        stack = []
        pairs = []
        print(dot_bracket)
        
        for i, char in enumerate(dot_bracket):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    pairs.append((stack.pop(), i))

        print("PAIRS FOUND BY PARSER:")
        print(pairs)
        
        return pairs

    # Function to calculate pair frequencies from a list of dot-bracket structures
    def calculate_pair_frequencies(dot_bracket_list):
        pair_counts = defaultdict(int)
        total_pairs = 0
        
        for dot_bracket in dot_bracket_list:
            pairs = parse_dot_bracket(dot_bracket)
            for pair in pairs:
                pair_counts[pair] += 1
                total_pairs += 1
        
        pair_frequencies = {pair: count / total_pairs for pair, count in pair_counts.items()}
        print("PAIR FREQUENCIES:")
        print(pair_frequencies)
        
        return pair_frequencies

    # Function to map nucleotide to color
    def nucleotide_color(nucleotide):
        color_map = {
            'A': 'red',
            'U': 'green',
            'G': 'blue',
            'C': 'orange'
        }
        return color_map.get(nucleotide, 'black')  # Default to black if unknown nucleotide

    # Function to plot arc diagram
    def plot_arc_diagram(sequence, pair_frequencies, pairs):
        
        print(pairs)
        
        fig = go.Figure()
        
        # Plot nodes with nucleotide-specific colors
        for i, nucleotide in enumerate(sequence):
            color = nucleotide_color(nucleotide)
            fig.add_trace(go.Scatter(
                x=[i], y=[0], mode='markers+text', text=nucleotide, textposition="bottom center",
                marker=dict(size=10, color=color),
                showlegend=False
            ))
        
        # Get the frequency for each pair
        pairs_unique = list(pair_frequencies.keys())
        frequencies = [pair_frequencies.get(pair, 0) for pair in pairs_unique]

        min_frequency = min(frequencies)
        max_frequency = max(frequencies)
        
        if frequencies:
            norm_frequencies = [(f - min(frequencies)) / (max(frequencies) - min(frequencies)) for f in frequencies]
        else:
            norm_frequencies = [0] * len(pairs)
        
        # Custom color scale from blue (least frequent) to red (most frequent)
        custom_colorscale = [
            [0.0, 'rgb(0, 0, 255)'],   # Blue
            [1.0, 'rgb(255, 0, 0)']    # Red
        ]
        
        # Plot arcs for edges with frequency-based colors
        for (start, end), norm_frequency, frequency in zip(pairs_unique, norm_frequencies, frequencies):
            center = (start + end) / 2
            radius = (end - start) / 2
            theta = np.linspace(0, np.pi, 100)
            x = center + radius * np.cos(theta)
            y = radius * np.sin(theta)
            color = px.colors.sample_colorscale(custom_colorscale, norm_frequency)
            fig.add_trace(go.Scatter(
                x=x, y=y, mode='lines',
                line=dict(color=color[0], width=2),
                showlegend=False,
                hoverinfo='text',
                text=f'Frequency: {frequency:.2f}'
            ))
        
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(
                colorscale=custom_colorscale,
                cmin=min_frequency, cmax=max_frequency,
                colorbar=dict(
                    title="Frequency",
                    tickvals=np.linspace(min_frequency, max_frequency, num=5),
                    ticktext=[f'{v:.2f}' for v in np.linspace(min_frequency, max_frequency, num=5)],
                    x=1.1,  # Position color bar to the right of the plot
                    y=0.5,
                    len=0.75,
                    thickness=20
                ),
                showscale=True
            ),
            hoverinfo='none'
        ))
        
        # Update layout
        fig.update_layout(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            showlegend=False,
            width=800,
            height=400,
            plot_bgcolor='white',
            title="RNA Secondary Structure Arc Diagram with Frequency-Based Coloring"
        )
        return fig
    
    print(dot_bracket_list)
    # Calculate pair frequencies from the dot-bracket list
    pair_frequencies = calculate_pair_frequencies(dot_bracket_list)
    print("FINAL PAIR FREQUENCIES:")
    print(pair_frequencies)


    pairs = parse_dot_bracket(dotbracket_native[0])
    print("PAIRS FOUND NATIVE")
    print(pairs)
    # Plot arc diagram
    fig = plot_arc_diagram(sequence, pair_frequencies, pairs)
    return fig



@socketio.on("connection")
def handle_message(data):
    print(data)


@socketio.on("connect")
def handle_connect():
    session_id = request.args.get("session_id")
    if session_id:
        join_room(session_id)
        socketio.emit(
            "command_output",
            {"output": "Connected to room: " + session_id},
            room=session_id,
        )
    print("Client connected")


@socketio.on("disconnect")
def handle_disconnect():
    session_id = request.args.get("session_id")
    if session_id:
        leave_room(session_id)
    print("Client disconnected")

@app.route("/upload-files", methods=["POST"])
def upload_files():
    if "nativePdb" not in request.files or "trajXtc" not in request.files:
        return redirect(request.url)

    native_pdb = request.files["nativePdb"]
    traj_xtc = request.files["trajXtc"]

    if native_pdb.filename == "" or traj_xtc.filename == "":
        return redirect(request.url)

    app.logger.info(request.form)  # Log request.form to check for data

    session_id = uuid.uuid4().hex
    socketio.emit("command_output", {"output": "BONJOUR"}, room=session_id)

    os.makedirs(os.path.join(app.static_folder, session_id), exist_ok=True)

    native_pdb_path = os.path.join(
        app.static_folder, session_id, secure_filename(native_pdb.filename)
    )
    traj_xtc_path = os.path.join(
        app.static_folder, session_id, secure_filename(traj_xtc.filename)
    )

    native_pdb.save(native_pdb_path)
    traj_xtc.save(traj_xtc_path)
    app.logger.info("SAVED TRAJECTORY")

    selected_plots = []
    n_frames = int(request.form.get("n_frames", 1))  # Default to 1 if not specified
    frame_range = request.form.get(
        "frame_range", "all"
    )  # Default to 'all' if not specified

    for plot in [
        "RMSD",
        "ERMSD",
        "CONTACT_MAPS",
        "TORSION",
        "SEC_STRUCTURE",
        "DOTBRACKET",
        "ARC",
        "ANNOTATE",
        "DS_MOTIF",
        "SS_MOTIF",
        "JCOUPLING",
        "ESCORE",
        "LANDSCAPE",
    ]:  # Adjust based on your available plots
        app.logger.info(plot)
        print(str(plot.lower()))
        if plot.lower() in request.form:
            app.logger.info(request.form)
            app.logger.info("In the Request part of the code")
            print(plot)
            selected_plots.append(plot)

    session["selected_plots"] = selected_plots
    session["n_frames"] = n_frames
    session["frame_range"] = frame_range
    session["torsionResidue"] = request.form.get("torsionResidue", 0)
    session["landscape_stride"] = request.form.get("landscape_stride", 0)
    session["form"] = list(request.form)
    session["files"] = {
        "nativePdb": str(request.files["nativePdb"].filename),
        "trajXtc": str(request.files["trajXtc"].filename),
        "trajpdb": str(request.files["trajXtc"].filename)
    }

    app.logger.info(selected_plots)

    return redirect(
        url_for(
            "view_trajectory",
            session_id=session_id,
            native_pdb=native_pdb.filename,
            traj_xtc=traj_xtc.filename,
        )
    )


@app.route("/retrieve-results", methods=["GET"])
def retrieve_results():
    print("RETRIVE RESULTS")
    session_id = request.args.get("session_id")
    print(session_id)

    directory_path = os.path.join(app.static_folder, session_id)
    pickle_file_path = os.path.join(directory_path, "plot_data.pkl")

    with open(os.path.join(directory_path, "files.json"), "r") as file:
        files_data = json.load(file)
        print(files_data)

    native_pdb = files_data['nativePdb']
    traj_xtc = files_data['trajXtc']

    
    native_pdb_path = os.path.join(directory_path, native_pdb)
    traj_xtc_path = os.path.join(directory_path, traj_xtc)
    u = mda.Universe(native_pdb_path, traj_xtc_path)

    if not os.path.exists(pickle_file_path):
        return "Session results not found", 404

    with open(pickle_file_path, "rb") as f:
        plot_data = pickle.load(f)

    with open(os.path.join(app.static_folder,'explanations.json')) as f:
        explanations = json.load(f)
    
    print(explanations)

    return render_template(
        "view_trajectory.html",
        session_id=session_id,
        native_pdb=native_pdb,
        traj_xtc=traj_xtc,
        plot_data=plot_data,
        trajectory_length=len(u.trajectory),
        explainations=explanations
    )


@app.route("/download/plot_data/<session_id>/<plot_id>")
def download_plot_data(plot_id, session_id):
    directory_path = os.path.join(app.static_folder, session_id)
    download_path = os.path.join(directory_path, "download_data")
    files_path = os.path.join(download_path, plot_id)
    print(f"Download path: {files_path}")
    file_paths = glob.glob(f"{files_path}/*")
    print(file_paths)

    if len(file_paths) > 1:
        memory_file = io.BytesIO()
        download_filename = f"{plot_id}_data.zip"
        with zipfile.ZipFile(memory_file, "w", zipfile.ZIP_DEFLATED) as zf:
            for file_path in file_paths:
                with open(file_path, "rb") as f:
                    # Use the file name as the arcname
                    zf.writestr(file_path.split("/")[-1], f.read())
    else:
        memory_file = io.BytesIO()
        with open(file_paths[0], "rb") as f:
            memory_file.write(f.read())
        download_filename = os.path.basename(file_paths[0])

    memory_file.seek(0)

    # Send the zip file as a response
    return send_file(memory_file, download_name=download_filename, as_attachment=True)


@app.route("/download/plot/<session_id>/<plot_id>")
def download_plot(plot_id, session_id):
    directory_path = os.path.join(app.static_folder, session_id)
    download_path = os.path.join(directory_path, "download_plot")
    files_path = os.path.join(download_path, plot_id)
    print(f"Download path: {files_path}")
    file_paths = glob.glob(f"{files_path}/*")
    print(file_paths)

    if len(file_paths) > 1:
        memory_file = io.BytesIO()
        download_filename = f"{plot_id}_plots.zip"
        with zipfile.ZipFile(memory_file, "w", zipfile.ZIP_DEFLATED) as zf:
            for file_path in file_paths:
                with open(file_path, "rb") as f:
                    # Use the file name as the arcname
                    zf.writestr(file_path.split("/")[-1], f.read())
    else:
        memory_file = io.BytesIO()
        with open(file_paths[0], "rb") as f:
            memory_file.write(f.read())
        download_filename = os.path.basename(file_paths[0])

    memory_file.seek(0)

    # Send the zip file as a response
    return send_file(memory_file, download_name=download_filename, as_attachment=True)


@app.route("/view-trajectory/<session_id>/<native_pdb>/<traj_xtc>")
def view_trajectory(session_id, native_pdb, traj_xtc):
    print("LOADED VIEW_TRAJECTORY")
    socketio.emit("command_output", {"output": "BONJOUR"}, room=session_id)
    directory_path = os.path.join(app.static_folder, session_id)
    download_path = os.path.join(directory_path, "download_data")
    download_plot = os.path.join(directory_path, "download_plot")

    with open(os.path.join(app.static_folder,'explanations.json')) as f:
        explanations = json.load(f)

    if not os.path.exists(download_path):
        os.makedirs(download_path)
    if not os.path.exists(download_plot):
        os.makedirs(download_plot)

    if not os.path.exists(directory_path):
        return "Session not found", 404
    
    print(f"FORM: {session['files']}")

    with open(os.path.join(directory_path, "files.json"), "w") as json_file:
        json.dump(session['files'], json_file)

    native_pdb_path = os.path.join(directory_path, native_pdb)
    traj_xtc_path = os.path.join(directory_path, traj_xtc)
    if not os.path.exists(native_pdb_path) or not os.path.exists(traj_xtc_path):
        return "File not found", 404

    selected_plots = session.get("selected_plots", [])
    n_frames = session.get("n_frames", 1)
    frame_range = session.get("frame_range", "all")
    app.logger.info(selected_plots)

    # Load trajectory with frame selection
    u = mda.Universe(native_pdb_path, traj_xtc_path)
    if frame_range != "all":
        start, end = map(int, frame_range.split("-"))
        u.trajectory[start:end:n_frames]
    else:
        u.trajectory[::n_frames]

    plot_data = []
    print(f"TRAJECTORY: {len(u.trajectory)}")
    
    
    image = generate_image(100)
    
    # Save the image to a BytesIO object
    img_io = BytesIO()
    image.save(img_io, 'PNG')
    img_io.seek(0)
    
    # Save the image to the static folder
    image.save('static/generated_image.png')


    socketio.emit("command_output", {"output": "LOUIS_TEST"})
    socketio.emit("command_output", {"output": selected_plots[0]})
    for plot in selected_plots:
        app.logger.info(plot)
        files_path = os.path.join(download_path, plot)
        plot_path = os.path.join(download_plot, plot)
        if not os.path.exists(files_path):
            os.makedirs(files_path)

        if not os.path.exists(plot_path):
            os.makedirs(plot_path)

        if plot in ["DOTBRACKET", "ARC"]:
            stackings, pairings, res = bb.annotate(
                traj_xtc_path, topology=native_pdb_path
            )
            dotbracket_data = bb.dot_bracket(pairings, res)[0]

        if plot == "CONTACT_MAPS":
            create_contact_maps(traj_xtc, native_pdb, 1, session_id)
            plot_data.append([plot, "CONTACT_MAPS", "rien"])


        if plot == "RMSD":
            rmsd = bb.rmsd(
                native_pdb_path,
                traj_xtc_path,
                topology=native_pdb_path,
                heavy_atom=True,
            )
            fig = plot_rmsd(rmsd)
            rmsd_df = pd.DataFrame({"RMSD": rmsd})
            rmsd_df.to_csv(os.path.join(files_path, "rmsd_values.csv"), index=False)
            fig.write_html(os.path.join(plot_path, "rmsd_plot.html"))

            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "scatter2", plotly_data])
        elif plot == "ERMSD":
            ermsd = bb.ermsd(native_pdb_path, traj_xtc_path, topology=native_pdb_path)
            fig = plot_ermsd(ermsd)
            ermsd_df = pd.DataFrame({"ERMSD": ermsd})
            ermsd_df.to_csv(os.path.join(files_path, "ermsd_values.csv"), index=False)
            fig.write_html(os.path.join(plot_path, "ermsd_plot.html"))

            plotly_data = plotly_to_json(fig)
            
            plot_data.append([plot, "scatter", plotly_data])
            app.logger.info(plotly_data)
        elif plot == "TORSION":
            angles, res = bb.backbone_angles(traj_xtc_path, topology=native_pdb_path)
            print(res)
            fig = plot_torsion(angles, res, session["torsionResidue"])
            fig.write_html(os.path.join(plot_path, "torsion_plot.html"))
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "torsion", plotly_data])
        elif plot == "SEC_STRUCTURE":
            sec_structure = bb.secondary_structure(
                traj_xtc_path, topology=native_pdb_path
            )
            fig = plot_sec_structure(sec_structure)
            fig.write_html(os.path.join(plot_path, "secondary_structure_plot.html"))
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "sec_structure", plotly_data])
        elif plot == "ANNOTATE":
            stackings, pairings, res = bb.annotate(
                traj_xtc_path, topology=native_pdb_path
            )
            plot_data.append([plot, "annotate", stackings, pairings, res])
        elif plot == "DOTBRACKET":
            fig = plot_dotbracket(dotbracket_data)
            fig.write_html(os.path.join(plot_path, "dotbracket_timeline_plot.html"))
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "dotbracket", plotly_data])
        elif plot == "ARC":
            print("ARC")

            #Dotbracket for the native state:
            stackings, pairings, res = bb.annotate(native_pdb_path)
            dotbracket_native = bb.dot_bracket(pairings, res)[0]
            print(f"DOTBRACKET NATIVE = {str(dotbracket_native[0])}")

            sequence = ''.join([item[0] for item in res])
            dotbracket_df = pd.DataFrame(dotbracket_data, columns=["DotBracket"])
            dotbracket_df.to_csv(os.path.join(files_path, "dotbracket_data.csv"), index=False)
            fig = plot_diagram_frequency(sequence, dotbracket_data, dotbracket_native)
            fig.write_html(os.path.join(plot_path, "arc_diagram_plot.html"))
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "arc", plotly_data])



        elif plot == "DS_MOTIF":
            # Implement DS_MOTIF plot
            pass
        elif plot == "SS_MOTIF":
            # Implement SS_MOTIF plot
            pass
        elif plot == "JCOUPLING":
            couplings, res = bb.jcouplings(traj_xtc_path, topology=native_pdb_path)
            plot_data.append([plot, "jcoupling", couplings])
        elif plot == "ESCORE":
            # Implement ESCORE plot
            pass

        elif plot == "LANDSCAPE":
            trajectory = Trajectory(
                filename=traj_xtc_path, ref_filename=native_pdb_path
            )
            test = trajectory.q_soft
            values = []
            good_rmsd = []
            rmsd = bb.rmsd(
                native_pdb_path,
                traj_xtc_path,
                topology=native_pdb_path,
                heavy_atom=True,
            )
            print(len(rmsd))

            for x in range(len(trajectory)):
                if x % int(session["landscape_stride"]) == 0:
                    print(x)
                    values.append(test[x])
                    good_rmsd.append(rmsd[x])

            df = pd.DataFrame(
                {
                    "frame": list(range(0, len(good_rmsd))),
                    "Q": values,
                    "RMSD": good_rmsd,
                    "traj": "traj_1",
                }
            )
            print("finished dataframe")
            size = 65
            # selected_regions=[(0.06,0.17,0.95,1.05),(0.15,0.22,0.78,0.92),(0.2,0.27,0.67,0.77),(0.25,0.33,0.49,0.63)] #list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!
            selected_regions = (
                []
            )  # list of tuples with minQ, maxQ, minRMSD, maxRMSD FROM CLOSEST TO NATIVE TO UNFOLDED!!!

            dataframe, max_RMSD = df, max(df["RMSD"])
            (
                probability_matrix,
                allframes_matrix,
                Qbin,
                RMSDbin,
            ) = energy_3dplot.make_matrix_probability(dataframe, size, max_RMSD)
            # energy_3dplot.probability_plot(probability_matrix, Qbin, RMSDbin, max_RMSD)
            energy_matrix, real_values = energy_3dplot.make_matrix_energy(
                probability_matrix, max_RMSD, size
            )
            fig = plot_landscapes_3D(
                energy_matrix, Qbin, RMSDbin, max_RMSD, real_values, selected_regions
            )

            path_landscape = f"static/{session_id}/landscape.png"
            energy_3dplot.energy_plot(
                energy_matrix,
                Qbin,
                RMSDbin,
                max_RMSD,
                real_values,
                selected_regions,
                path_landscape,
            )
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "landscape", plotly_data])
            path_landscape = f"{session_id}/landscape.png"
            plot_data.append(["LANDSCAPE_PLT", "matplotlib", path_landscape])

        pickle_file_path = os.path.join(directory_path, "plot_data.pkl")
        with open(pickle_file_path, "wb") as f:
            pickle.dump(plot_data, f)

    return render_template(
        "view_trajectory.html",
        session_id=session_id,
        native_pdb=native_pdb,
        traj_xtc=traj_xtc,
        plot_data=plot_data,
        trajectory_length=len(u.trajectory),
        explainations=explanations,
    )


@app.route("/plot-biased", methods=["POST"])
def plot_biased():
    directory = request.form["directory"]
    df = read_biased(os.path.join(directory, "ratchet.out"))
    plot_filenames = save_all_plots(df, directory)
    return render_template("biased_plots.html", plot_filenames=plot_filenames)


def read_biased(file_path):
    df = pd.read_csv(file_path, sep="\t", header=None)
    return df

def save_nth_frame(traj_file, top_file, index, session):

    # Load the trajectory and topology
    u = mda.Universe(top_file, traj_file)

    # Select the nth frame
    u.trajectory[index]

    # Define the path to save the frame
    path_save = os.path.join('static', session, 'temp_frame.pdb')

    # Save the frame
    with mda.Writer(path_save, multiframe=False) as W:
        W.write(u)

    return path_save

def create_contact_maps(traj_file, top_file, index, session):
    path_save_contact_maps = os.path.join('static', session, 'target.cmp')
    path_save_figure = os.path.join('static', session, 'contact_map_plotly.png')

    path_traj_file = os.path.join('static', session, traj_file)
    path_top_file = os.path.join('static', session, top_file)

    path_save_nth_frame = save_nth_frame(path_traj_file, path_top_file, index, session)
    contact_map, signature = generate_contact_maps.get_contact_map_and_signature(path_save_nth_frame)
    generate_contact_maps.write_contact_map(path_save_contact_maps, contact_map, signature)
    generate_contact_maps.read_contact_map(path_save_contact_maps, path_save_figure)




@socketio.on('slider_value')
def handle_slider_value(data):
    slider_value = int(data['value'])
    traj = data['traj']
    native_state = data['native_state']
    session = data['session_id']
    print(f"slider value = {slider_value}")
    path_contactmap_figure = create_contact_maps(traj, native_state, slider_value, session)


    path_save_2 = 'contact_map_plotly.png'

    emit('image_update', {'image_url': path_save_2})

def generate_image(value):
    # Create a simple image based on the slider value
    print("generating image")
    width, height = 200, 100
    image = Image.new('RGB', (width, height), color=(255, 255, 255))
    draw = ImageDraw.Draw(image)
    draw.rectangle([10, 10, 10 + value, 90], fill=(255, 0, 0))
    return image


def save_plot(df, time_column, column_to_plot, output_folder):
    plt.figure(figsize=(10, 6))
    plt.plot(df[time_column], df[column_to_plot])
    plt.title(f"{column_to_plot} over time")
    plt.xlabel("Time")
    plt.ylabel(column_to_plot)
    plot_filename = f"{column_to_plot}_over_time.png"
    plt.savefig(os.path.join(output_folder, plot_filename))
    plt.close()
    return plot_filename


def save_all_plots(df, directory):
    output_folder = os.path.join(directory, "plots")
    os.makedirs(output_folder, exist_ok=True)
    time_column = df.columns[0]
    plot_filenames = []
    for column_to_plot in df.columns[1:]:
        plot_filename = save_plot(df, time_column, column_to_plot, output_folder)
        plot_filenames.append(plot_filename)
    return plot_filenames


if __name__ == "__main__":
    # app.run(debug=True)
    socketio.run(app)
