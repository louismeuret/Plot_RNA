from flask import request, render_template, redirect, url_for, send_file
from flask import Flask, request, render_template, send_from_directory
import plotly.graph_objs as go

app = Flask(__name__)
# Dictionary containing plot options
plot_options = {
    'plot1': 'Plot 1',
    'plot2': 'Plot 2',
    'plot3': 'Plot 3',
    'plot4': 'Plot 4',
    'plot5': 'Plot 5',
    'plot6': 'Plot 6',
    'plot7': 'Plot 7',
    'plot8': 'Plot 8'
}

@app.route('/generate-plots', methods=['POST'])
def generate_plots():
    # Retrieve form data
    structure1 = request.form.get('structure1')
    structure2 = request.form.get('structure2')

    # Get selected plot options
    selected_plots = [key for key in plot_options.keys() if request.form.get(key) == key]

    # Generate plots
    plots = {}
    for plot_id in selected_plots:
        if plot_id == 'plot1':
            # Example plot 1: scatter plot
            x_data = [1, 2, 3, 4, 5]
            y_data = [2, 3, 4, 5, 6]
            fig = go.Figure(go.Scatter(x=x_data, y=y_data, mode='lines+markers', name='Plot 1'))
        elif plot_id == 'plot2':
            # Example plot 2: bar chart
            categories = ['A', 'B', 'C', 'D', 'E']
            values = [10, 20, 15, 25, 30]
            fig = go.Figure(go.Bar(x=categories, y=values, name='Plot 2'))
        elif plot_id == 'plot3':
            # Add more plot options as needed
            pass
        # Add more conditions for other plots

        # Convert plot to HTML string and store
        plots[plot_id] = fig.to_html(full_html=False)

    # Store plots in session or pass to the template
    session['plots'] = plots

    # Redirect back to index
    return redirect(url_for('index'))

@app.route('/download-results')
def download_results():
    # Generate downloadable HTML file containing the plots
    plots = session.get('plots', {})
    if not plots:
        # No plots generated
        return "No plots generated yet."

    # Combine plots into a single HTML file
    combined_html = '<html><head><title>Plots</title></head><body>'
    for plot_name, plot_html in plots.items():
        combined_html += f'<h2>{plot_name}</h2>{plot_html}'
    combined_html += '</body></html>'

    # Save combined HTML to a temporary file
    temp_file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'plots.html')
    with open(temp_file_path, 'w') as f:
        f.write(combined_html)

    # Send the file to the user for download
    return send_file(temp_file_path, as_attachment=True)
