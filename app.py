from flask import Flask, request, render_template, send_from_directory,redirect,url_for, session
from werkzeug.utils import secure_filename
import uuid
import os
import glob
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
from MDAnalysis.topology.guessers import guess_types
from MDAnalysis.analysis.align import alignto
import json
import barnaba as bb
import plotly
from plotly.tools import mpl_to_plotly

app = Flask(__name__)
app.secret_key = "pi"

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        directory = request.form['directory']
        # Perform actions based on the form input, e.g., list files
        files = list_files(directory)
        return render_template('results.html', files=files)
    return render_template('index.html')

def list_files(directory):
    # Your implementation of listing files
    return os.listdir(directory)  # Simplified example

def plotly_to_json(fig):
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

def plot_ermsd(ermsd):
    fig = go.Figure(data=go.Scatter(y=ermsd, mode='markers'))
    fig.update_layout(title='ERMSD Scatter Plot', xaxis_title='Frame', yaxis_title='ERMSD Value')
    return fig

def plot_rmsd(rmsd):
    fig = go.Figure(data=go.Scatter(y=rmsd, mode='markers'))
    fig.update_layout(title='RMSD Scatter Plot', xaxis_title='Frame', yaxis_title='RMSD Value')
    return fig

def plot_dotbracket(dotbracket_data):
    reverse_mapping = {}
    for frame, dotbracket in enumerate(dotbracket_data, start=1):
        reverse_mapping.setdefault(dotbracket, []).append(frame)

    # Flatten the data structure and include frame information
    structures_data = []
    for dotbracket, frames in reverse_mapping.items():
        for frame in frames:
            structures_data.append({
                'Frame': frame,
                'DotBracket': dotbracket
            })

    # Create a DataFrame from the flattened data
    structures_df = pd.DataFrame(structures_data)

    # Sort the DataFrame by structure and then by frame
    structures_df.sort_values(by=['DotBracket', 'Frame'], inplace=True)

    # Get unique dot-bracket structures for color coding
    unique_structures = structures_df['DotBracket'].unique()
    structure_colors = [f'rgb{tuple(int(255 * x) for x in plt.cm.tab20b(i)[:3])}' for i in np.linspace(0, 1, len(unique_structures))]
    color_dict = dict(zip(unique_structures, structure_colors))

    # Create traces for each structure
    traces = []
    for i, dotbracket in enumerate(unique_structures):
        structure_df = structures_df[structures_df['DotBracket'] == dotbracket]
        x_values = []
        y_values = []
        # Track the previous frame for comparison
        prev_frame = None
        for _, row in structure_df.iterrows():
            frame = row['Frame']
            if prev_frame is not None and frame != prev_frame + 1:
                # Add None to create gaps between structures
                x_values.append(None)
                y_values.append(None)
            x_values.append(frame)
            y_values.append(i + 1 - 0.2)  # Adjust the y-value slightly for better visualization
            prev_frame = frame
        trace = go.Scatter(
            x=x_values,
            y=y_values,
            mode='lines',
            line=dict(color=color_dict[dotbracket], width=8),
            name=dotbracket
        )
        traces.append(trace)

    # Layout
    layout = go.Layout(
        title='Timeline of RNA Structures',
        xaxis=dict(
            title='Frame',
            showgrid=False,
            zeroline=False,
            showline=False
        ),
        yaxis=dict(
            title='Dot-Bracket Structure',
            tickmode='array',
            tickvals=list(range(1, len(unique_structures) + 1)),
            ticktext=unique_structures,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticklen=1
        ),
        showlegend=False,
        plot_bgcolor='rgba(255, 255, 255, 0)'
    )

    # Plot
    fig = go.Figure(data=traces, layout=layout)
    return fig

def plot_torsion(torsion):
    return fig 

def plot_sec_structure(sec_structure):
    return fig 







@app.route('/upload-files', methods=['POST'])
def upload_files():
    if 'nativePdb' not in request.files or 'trajXtc' not in request.files:
        return redirect(request.url)
    
    native_pdb = request.files['nativePdb']
    traj_xtc = request.files['trajXtc']
    
    if native_pdb.filename == '' or traj_xtc.filename == '':
        return redirect(request.url)
    
    app.logger.info(request.form)  # Log request.form to check for data
    
    # Generate a unique session ID for this user's files
    session_id = uuid.uuid4().hex
    
    os.makedirs(os.path.join(app.static_folder, session_id), exist_ok=True)
    
    native_pdb_path = os.path.join(app.static_folder, session_id, secure_filename(native_pdb.filename))
    traj_xtc_path = os.path.join(app.static_folder, session_id, secure_filename(traj_xtc.filename))
    
    native_pdb.save(native_pdb_path)
    traj_xtc.save(traj_xtc_path)
    app.logger.info("SAVED TRAJECTORY")
    
    selected_plots = []

    for plot in ['RMSD','ERMSD', 'TORSION', 'SEC_STRUCTURE','DOTBRACKET', 'ANNOTATE', 'DS_MOTIF', 'SS_MOTIF','JCOUPLING','ESCORE']:  # Adjust based on your available plots
        app.logger.info(plot)
        if plot.lower() in request.form:
            app.logger.info(request.form)
            app.logger.info("In the Request part of the code")
            selected_plots.append(plot)
                
    
    session['selected_plots'] = selected_plots
    app.logger.info(selected_plots)
    return redirect(url_for('view_trajectory', session_id=session_id, native_pdb=native_pdb.filename, traj_xtc=traj_xtc.filename))

@app.route('/view-trajectory/<session_id>/<native_pdb>/<traj_xtc>')
def view_trajectory(session_id, native_pdb, traj_xtc):
    # Ensure the directory and files exist
    directory_path = os.path.join(app.static_folder, session_id)
    if not os.path.exists(directory_path):
        return "Session not found", 404
    
    native_pdb_path = os.path.join(directory_path, native_pdb)
    traj_xtc_path = os.path.join(directory_path, traj_xtc)
    if not os.path.exists(native_pdb_path) or not os.path.exists(traj_xtc_path):
        return "File not found", 404
    
    # Retrieve selected plots from session data
    selected_plots = session.get('selected_plots', [])
    app.logger.info(selected_plots)
    plot_data = []
    for plot in selected_plots:  # Adjust based on your available plots
        print(plot)
        app.logger.info(plot)
        if plot == 'RMSD':
            # Perform specific action for RMSD plot
            #stackings, pairings, res = bb.annotate(traj_xtc.filename,native_pdb.filename)
            rmsd = bb.rmsd(native_pdb_path,traj_xtc_path,topology=native_pdb_path, heavy_atom=True)
            fig = plot_rmsd(rmsd)
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "scatter2", plotly_data])

        elif plot == 'ERMSD':
            # Need to be no accent in the path, else it crash
            ermsd = bb.ermsd(native_pdb_path,traj_xtc_path,topology=native_pdb_path)
            fig = plot_ermsd(ermsd)
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "scatter", plotly_data])
            app.logger.info(plotly_data)

        elif plot == 'TORSION':
            # Perform specific action for TORSION plot
            angles, res = bb.backbone_angles(traj, topology=top)
            selected_plots.append(plot)
        elif plot == 'SEC_STRUCTURE':
            # Perform specific action for SEC_STRUCTURE plot
            selected_plots.append(plot)
        elif plot == 'ANNOTATE':
            # Perform specific action for ANNOTATE plot
            stackings, pairings, res = bb.annotate(traj, topology=top)
            selected_plots.append(plot)
        elif plot == 'DOTBRACKET':
            stackings, pairings, res = bb.annotate(traj_xtc_path,topology=native_pdb_path)
            dotbracket_data = bb.dot_bracket(pairings,res)[0]
            fig = plot_dotbracket(dotbracket_data)
            plotly_data = plotly_to_json(fig)
            plot_data.append([plot, "dotbracket", plotly_data])

        elif plot == 'DS_MOTIF':
            selected_plots.append(plot)
        elif plot == 'SS_MOTIF':
            # Perform specific action for SS_MOTIF plot
            selected_plots.append(plot)
        elif plot == 'JCOUPLING':
            # Perform specific action for JCOUPLING plot
            couplings, res = bb.jcouplings(traj, topology=top)
            selected_plots.append(plot)
        elif plot == 'ESCORE':
            # Perform specific action for ESCORE plot
            selected_plots.append(plot)
        elif plot == 'DOT_BRACKET':
            dotbr, ss = bb.dot_bracket(pairings, res)
        elif plot == 'GVEC':
            gvec, seq = bb.dump_gvec(traj, topology=top)

    print(plotly_data)


    
    # Generate or fetch plot filenames
    plot_filenames = [f'plot_{plot}.png' for plot in selected_plots]
    return render_template('view_trajectory.html', session_id=session_id, native_pdb=native_pdb, traj_xtc=traj_xtc, plot_filenames=plot_filenames,graphJSON=plotly_data, plot_data=plot_data)



@app.route('/plot-biased', methods=['POST'])
def plot_biased():
    directory = request.form['directory']
    df = read_biased(os.path.join(directory, "ratchet.out"))
    # Assuming a function that saves a plot and returns its filename
    plot_filename = save_plot(df)
    return send_from_directory(directory='static', filename=plot_filename)

def read_biased(file_path):
    df = pd.read_csv(file_path, sep="\t", header=None)
    # Perform your data manipulation
    return df

def save_plot(df):
    # Plotting logic, save plot to static folder, and return its filename
    plt.figure()
    df.plot()
    plot_filename = 'plot.png'
    plt.savefig(f'static/{plot_filename}')
    plt.close()
    return plot_filename



# Additional routes for other functionalities as needed

if __name__ == '__main__':
    app.run(debug=True)

