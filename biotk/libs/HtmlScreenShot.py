import os
from selenium import webdriver

def which(exe):
    """
    Find path to specified executable
    """
    PATH = os.getenv('PATH')
    for place in PATH.split(os.path.pathsep):
        place = os.path.join(place, exe)
        if os.path.exists(place) and os.access(place, os.X_OK):
            return place
    return None

class PhantomDriver:
    """
    Class for initializing a PhantomJS executable driver
    """
    def __init__(self, exe):
        self.exe = exe
        self.phantomjs_river = None

    def __enter__(self):
        self.phantomjs_driver = webdriver.PhantomJS(executable_path=self.exe)
        return self.phantomjs_driver

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.phantomjs_driver is not None:
            self.phantomjs_driver.quit()

class HtmlScreenShot:
    """
    Use a headless browser to take a screenshot of a
    webpage and save as static image.

    Initial usage for generating thumbnails of Plotly
    images for Zia
    """
    def __init__(self, html_path,
                       output_path=None):
        self.html_path = html_path
        self.output_path = output_path
        self.PH_EXE = which('phantomjs')
        if self.PH_EXE is not None:
            # take the screenshot
            ph_driver = PhantomDriver(self.PH_EXE)
            self.png_path = self.take_screenshot(ph_driver)
        else:
            raise OSError('Could not resolve path to PhantomJS executable')

    def take_screenshot(self, phantomjs_driver):
        """
        Use PhantomJS and Selenium to take screenshot
        of HTML page and save to PNG
        """
        if self.output_path is None:
            png_path = os.path.splitext(self.html_path)[0] + '.png'
        else:
            png_path = self.output_path
        with phantomjs_driver as driver:
            driver.set_window_size(1920, 1080)
            driver.get(self.html_path)
            driver.save_screenshot(png_path)

        return png_path
