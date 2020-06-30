import kivy
kivy.require('1.0.7')

from kivy.app import App

from kivy.uix.label import Label
from kivy.uix.gridlayout import GridLayout
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button

from kivy.garden.matplotlib.backend_kivyagg import FigureCanvasKivyAgg
from kivy.garden.filebrowser import FileBrowser
from kivy.uix.filechooser import FileChooserListView

from kivy.uix.popup import Popup

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import butter, lfilter
from astropy.table import Table

Browsertxt=r"G:\Stopped Flow Exp\2019.10.10"
plt.plot([1,1,1])

class MyFileChooser(FileChooserListView):
    
    def on_submit(*args, **kwargs):
        global browser
        browser=args[1][0]
        print(browser)
    
class MyGrid(GridLayout):
    
    def DirectoryGet(self, instance,*args):
        global Browsertxt
        layout = GridLayout(cols = 1)
        self.fclv = MyFileChooser(filters= [lambda folder, filename: not filename.endswith('.sys')])
        layout.add_widget(self.fclv)
        
        closeButton = Button(text = "Close the pop-up")
        layout.add_widget(closeButton) 
        popup = Popup(title ='Demo Popup',content = layout,auto_dismiss=True)   
        popup.open()
        closeButton.bind(on_press = popup.dismiss) 
        self.Directory=print(Browsertxt)
        return self.Directory
    def BufferGet(self, instance):
        print("pressed BufferGet")
    def QuencherGet(self, instance):
        print("pressed QuencherGet")
    def Standard(self, instance):
        Directory=self.Directory.text
        Buffer=self.Buffer.text
        BufferR=self.BufferR.text
        Quencher=self.Quencher.text
        QuencherR=self.QuencherR.text
        print("pressed Standard"+Directory+Buffer+Quencher+BufferR+QuencherR)
        plt.clf()
        plt.plot([0,1,2])
        self.remove_widget(self.Figure)
        self.Figure=FigureCanvasKivyAgg(plt.gcf())
        self.add_widget(self.Figure)
        
    def KSV(self, instance):
        Directory=self.Directory.text
        Buffer=self.Buffer.text
        BufferR=self.BufferR.text
        Quencher=self.Quencher.text
        QuencherR=self.QuencherR.text
        print("pressed KSV"+Directory+Buffer+Quencher+BufferR+QuencherR)
        plt.clf()
        plt.plot([2,1,0])
        self.remove_widget(self.Figure)
        self.Figure=FigureCanvasKivyAgg(plt.gcf())
        self.add_widget(self.Figure)
                
    def __init__(self, **kwargs):
        super(MyGrid,self).__init__(**kwargs)
        
        self.cols=1
        
        self.inside=GridLayout()
        self.inside.cols=4
        
        self.submit=Button(text="Directiory")
        self.submit.bind(on_press=self.DirectoryGet)
        self.inside.add_widget(self.submit)
        self.Directory=TextInput(text=Browsertxt,font_size=12)
        self.inside.add_widget(self.Directory)
        
        self.submit=Button(text="Standard Fit")
        self.submit.bind(on_press=self.Standard)
        self.inside.add_widget(self.submit)
        
        self.submit=Button(text="Stern Volmer Mod Fit")
        self.submit.bind(on_press=self.KSV)
        self.inside.add_widget(self.submit)
        
        self.submit=Button(text="Buffer")
        self.submit.bind(on_press=self.BufferGet)
        self.inside.add_widget(self.submit)
        self.Buffer=TextInput(text="Run00000.dsa")
        self.inside.add_widget(self.Buffer)
        
        self.Label=Label(text="Buffer Remover")
        self.inside.add_widget(self.Label)
        self.BufferR=TextInput(text="1 2")
        self.inside.add_widget(self.BufferR)
        
        self.submit=Button(text="Quencher")
        self.submit.bind(on_press=self.QuencherGet)
        self.inside.add_widget(self.submit)
        self.Quencher=TextInput(text="Run00001.dsa")
        self.inside.add_widget(self.Quencher)
        
        self.Label=Label(text="Quencher Remover")
        self.inside.add_widget(self.Label)
        self.QuencherR=TextInput(text="1 2")
        self.inside.add_widget(self.QuencherR)
        
        self.add_widget(self.inside)
        
        self.Figure=FigureCanvasKivyAgg(plt.gcf())
        self.add_widget(self.Figure)
        

        
class TestApp(App):
    def build(self):
        return MyGrid()


if __name__ == '__main__':
    TestApp().run()