#!/usr/bin/env python

import wx
class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent)

        #add position panel
        posPnl = wx.Panel(self)
        lbl1 = wx.StaticText(posPnl, label="Position")
        lbl2 = wx.StaticText(posPnl, label="Size")
        sizeCtrl = wx.TextCtrl(posPnl)

        posPnlSzr = wx.BoxSizer(wx.HORIZONTAL)
        posPnlSzr.Add(lbl1, 1, wx.GROW)
        posPnlSzr.Add(sizeCtrl, 1, wx.GROW)
        posPnlSzr.Add(lbl2, 1, wx.GROW)

        posPnl.SetSizer(posPnlSzr)

        #create a top leverl sizer to add to the frame itself
        mainSzr = wx.BoxSizer(wx.VERTICAL)
        mainSzr.Add(posPnl, 1, wx.GROW)

        self.SetSizerAndFit(mainSzr)
        self.Show()


app = wx.App(False)
frame = MainWindow(None, "Trading Client")
app.MainLoop()
