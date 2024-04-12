from flask import Flask, render_template, request
import plotly.express as px
import pandas as pd

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        f = request.files['file']
        df = pd.read_csv(f)
        fig = px.line(df, x="x_column", y="y_column")  # Adjust according to your data
        graphHTML = fig.to_html(full_html=False)
        return render_template("plot.html", plot=graphHTML)
    return render_template("upload.html")

if __name__ == '__main__':
    app.run(debug=True)

