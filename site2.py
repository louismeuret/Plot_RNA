from flask import Flask, render_template_string, jsonify
import random
import json

app = Flask(__name__)

@app.route('/')
def index():
    initial_data = {
        'x': [1, 2, 3, 4, 5],
        'y': [10, 11, 12, 13, 14],
        'type': 'scatter'
    }

    # Convert initial plot data to JSON
    initial_data_json = json.dumps(initial_data)

    # HTML template with added CSS for styling
    template = '''
    <!DOCTYPE html>
    <html>
        <head>
            <title>Interactive Plotly Plot</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                #plot { width: 100%; max-width: 600px; margin: auto; }
                button {
                    display: block;
                    margin: 20px auto;
                    padding: 10px 20px;
                    background-color: #007BFF;
                    border: none;
                    border-radius: 5px;
                    color: white;
                    cursor: pointer;
                    font-size: 16px;
                }
                button:hover {
                    background-color: #0056b3;
                }
                h1 { text-align: center; }
            </style>
        </head>
        <body onload="initPlot()">
            <h1>Interactive Plotly Plot in Flask</h1>
            <div id="plot"></div>
            <button onclick="requestData()">Request New Data</button>
            <script>
                const initialData = {{ initial_data_json|safe }};
                
                function initPlot() {
                    Plotly.newPlot('plot', [initialData], {responsive: true});
                }

                function requestData() {
                    fetch('/get-data')
                        .then(response => response.json())
                        .then(data => {
                            Plotly.react('plot', [data], {responsive: true});
                        });
                }
            </script>
        </body>
    </html>
    '''
    return render_template_string(template, initial_data_json=initial_data_json)

@app.route('/get-data', methods=['GET'])
def get_data():
    # Generate random data for the sake of example
    x = list(range(10))
    y = [random.randint(1, 20) for _ in x]
    plot_data = {
        'x': x,
        'y': y,
        'type': 'scatter'
    }
    return jsonify(plot_data)

if __name__ == '__main__':
    app.run(debug=True)
