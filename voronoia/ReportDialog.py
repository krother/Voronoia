
from Tkinter import *
import Pmw
import re,string,math

class ReportDialog(Pmw.Dialog):

    def __init__(self,parent,titles,labels,text,cmd):
        Pmw.Dialog.__init__(self,parent,title=titles,
                buttons = ['Save as HTML','Save as ASCII','OK'],
                command = cmd,
                )
        
        self.st = Pmw.ScrolledText(self.interior(),
                labelpos = 'n',
                label_text = labels,
                text_font = Pmw.logicalfont('Fixed'),
                text_wrap='none',

                usehullsize = 1,
                hull_width = 800,
                hull_height = 600,
                text_padx = 10,
                text_pady = 10,                
                
                )

        self.st.settext(text)
        self.st.pack(padx=10,pady=10)

        self.st.configure(
            text_state = 'disabled',
        )

    

