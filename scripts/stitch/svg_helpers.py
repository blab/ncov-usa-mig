from lxml import etree
from svgutils.compose import Element

class Rect(Element):
    """
    A custom SVG rectangle class for svgutils.compose.
    Use this to create backgrounds or colored blocks in stitched figures.
    """
    def __init__(self, width, height, fill="white"):
        self.root = etree.Element("rect", 
                                  width=f"{width}pt", 
                                  height=f"{height}pt", 
                                  fill=fill)