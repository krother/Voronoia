# -*- coding: cp1252 -*-

import re,string

class Report:
    """
    Class for converting HTML to text reports.
    """
    html_header = """
<html><head></head>
<body>
<p>   Report created by Voronoia 1.0<br>
   (c) Kristian Rother, Bjoern Gruening, Andrean Goede, Peter Hildebrand, Robert Preissner 2008<br>
   http://bioinf.charite.de
</p>
<hr>
    """
    html_footer = """
    <p><br></p>
<hr>
<p>
When using Voronoia for packing calculations (or anything else),
please cite:<br>
    Voronoia: Analyzing packing in protein structures<br>
    Kristian Rother, Bjoern Gruening, Andrean Goede, Peter Hildebrand, Robert Preissner (submitted)<br>
</p>
&nbsp;
</body></html>
    """

    def __init__(self):
        self.html = ""
        self.text = ""

    def add_html_body(self,html):
        self.html += html

    def get_html(self):
        return self.html_header+self.html+self.html_footer

    def get_text(self):
        self.convert()
        return self.text

    def convert(self):
        text = self.html+self.html_footer
        text = re.sub('\s*[\r\n]+','',text)
        text = re.sub('<table[^>]*>','\n',text)
        text = re.sub('<p>|</sup>|<i>|</i>|<b>|</b>|<td[^>]*>|<tr[^>]*>|<th[^>]*>|</li>','',text)
        text = re.sub('<h\d>|</body>|</html>|</table>|</tr>|<br>|</p>|<ul>|</ul>|<nl>|</nl>|&nbsp;','\n',text)
        text = re.sub('</h\d>','\n\n',text)
        text = re.sub('<hr>','----------------------------------------------------------------------\n',text)
        text = re.sub('<li>','\n * ',text)
        text = re.sub('</th>|</td>','\t',text)
        text = re.sub('&gt;','>',text)
        text = re.sub('&lt;','<',text)
        text = re.sub('<sup>','^',text)
        text = re.sub('&Aring;','Å',text)

        self.text = text
    
