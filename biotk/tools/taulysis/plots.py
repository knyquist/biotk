import os
import numpy as np
import plotly
from plotly.graph_objs import *
from plotly.offline import plot

def plot_line_plot(x,
                   y,
                   output_directory,
                   title,
                   datapoint_mode,
                   xlabel,
                   ylabel,
                   yaxis_range,
                   descriptors):
    """Use plotly to plot line or scatter plot"""
    data = Scatter(
        x=x,
        y=y,
        mode=datapoint_mode,
        marker=dict(
            opacity=0.25
        )
    )
    layout = Layout(
        title=title,
        font=dict(
            size=20
        ),
        xaxis=dict(
            title=xlabel#'Template position (bp)'
        ),
        yaxis=dict(
            title=ylabel,#'Fraction survivors left',
            range=yaxis_range
        )
    )
    figure = Figure(data=[data], layout=layout)
    additional_info = '-'.join([item for item in descriptors])
    filename = (additional_info + '_' + str(title) +
                '.html')
    filepath = os.path.sep.join([output_directory, filename])
    plot(figure, show_link=False, auto_open=False, filename=filepath)

def plot_fitted_line_plot(x,
                          y,
                          t_1,
                          t_2,
                          t_rc,
                          output_directory,
                          title,
                          xlabel,
                          ylabel,
                          descriptors):
    """Use plotly to plot line plot with horizontal fits"""
    data1 = Scatter(
        x=x,
        y=y,
        showlegend=False
    )
    # Fits
    # handle case where fit didn't happend
    if t_1[2] == 1:
        t_1[2] = 0
    if t_2[2] == 1:
        t_2[2] = 0
    if t_rc[2] == 1:
        t_rc[2] = 0

    data2 = Scatter(
        x=[t_1[0], t_1[1]],
        y=[t_1[2], t_1[2]],
        mode='lines',
        name=str(np.around(t_1[2], 7))
    )
    data3 = Scatter(
        x=[t_2[0], t_2[1]],
        y=[t_2[2], t_2[2]],
        mode='lines',
        name=str(np.around(t_2[2], 7))
    )
    data4 = Scatter(
        x=[t_rc[0], t_rc[1]],
        y=[t_rc[2], t_rc[2]],
        mode='lines',
        name=str(np.around(t_rc[2], 7))
    )
    data = [data1, data2, data3, data4]
    layout = Layout(
        title=title,
        font=dict(
            size=20
        ),
        xaxis=dict(
            title=xlabel
        ),
        yaxis=dict(
            title=ylabel
        )
    )
    figure = Figure(data=data, layout=layout)
    additional_info = '-'.join([item for item in descriptors])
    filename = (additional_info + '_' + str(title) +
                '.html')
    filepath = os.path.sep.join([output_directory, filename])
    plot(figure, show_link=False, auto_open=False, filename=filepath)

def plot_histogram(data,
                   output_directory,
                   title,
                   xlabel,
                   ylabel,
                   bins,
                   descriptors):
    """Use plotly to plot simple histogram"""
    data = Histogram(
        x=data,
        opacity=0.75,
        autobinx=False,
        xbins=dict(
            start=bins[0],
            end=bins[1],
            size=bins[2]
        )
    )
    layout = Layout(
        title=title,
        font=dict(
            size=20
        ),
        xaxis=dict(
            title=xlabel
        ),
        yaxis=dict(
            title=ylabel
        )
    )
    figure = Figure(data=[data], layout=layout)
    additional_info = '-'.join([item for item in descriptors])
    filename = (additional_info + '_' + str(title) +
                '.html')
    filepath = os.path.sep.join([output_directory, filename])
    plot(figure, show_link=False, auto_open=False, filename=filepath)
