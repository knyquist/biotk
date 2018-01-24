import biotk.libs.HtmlScreenShot as ss
import os

class _test_screen_shot:
    def __init__(self, html_path,
                       output_path):
        self.html_path = html_path
        self.output_path = output_path
        self.html_screen_shotter = ss.HtmlScreenShot(self.html_path,
                                                     self.output_path)

def test_screen_shot():
    """
    Given a Plotly HTML figure, demonstrate that
    a PNG file can be generated.
    """
    path = os.path.dirname(os.path.abspath(__file__))
    html_name = 'data/plotly_figure.html'
    html_path = os.path.join(path, html_name)
    png_name = 'plotly_figure.png'
    png_path = os.path.join('/tmp/', png_name)
    ss_test = _test_screen_shot(html_path,
                                png_path)

    assert ss_test.html_screen_shotter.png_path == png_path
