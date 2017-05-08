#!/usr/bin/python

import wx
import wx.lib.mixins.listctrl as listmix # for editable list ctrl

class EditableListCtrl(wx.ListCtrl, listmix.TextEditMixin):
	def __init__(self, parent, ID=wx.ID_ANY, pos=wx.DefaultPosition, size=wx.DefaultSize, style=0):
		wx.ListCtrl.__init__(self, parent, ID, pos, size, style)
		listmix.TextEditMixin.__init__(self)

class MainWindow(wx.Frame):
	def __init__(self, *args, **kwargs):
#	def __init__(self, parent, title):
		super(MainWindow, self).__init__(*args, **kwargs)
#		wx.Frame.__init__(self, parent, title=title, size=(200,100))

		sp = wx.SplitterWindow(self)

		spLeft = wx.SplitterWindow(sp)
		spRight = wx.SplitterWindow(sp)

		p1 = wx.Panel(spLeft, style=wx.SUNKEN_BORDER)
		p2 = wx.Panel(spLeft, style=wx.SUNKEN_BORDER)
		p3 = wx.Panel(spRight, style=wx.SUNKEN_BORDER)
		p4 = wx.Panel(spRight, style=wx.SUNKEN_BORDER)

		sp.SplitVertically(spLeft, spRight, 250)
		spLeft.SplitHorizontally(p1, p2, -200)
		spRight.SplitHorizontally(p3, p4, -200)

		sz1 = wx.BoxSizer(wx.HORIZONTAL)
		sz2 = wx.BoxSizer(wx.HORIZONTAL)
		sz3 = wx.BoxSizer(wx.HORIZONTAL)
		sz4 = wx.BoxSizer(wx.HORIZONTAL)

#		tp = TreePanel(p1)
#		tp = wx.TextCtrl(p1, 1, style=wx.TE_MULTILINE)
		tp = wx.TreeCtrl(p1, 1, wx.DefaultPosition, (-1,-1), wx.TR_HIDE_ROOT|wx.TR_HAS_BUTTONS)
		root = tp.AddRoot('Program')
		ga = tp.AppendItem(root, 'Genetic Algorithm')
		base = tp.AppendItem(ga, 'Base Parameters')
		mut = tp.AppendItem(ga, 'Mutation')
		b1 = tp.AppendItem(mut, 'Boolean')
		i1 = tp.AppendItem(mut, 'Integer')
		d1 = tp.AppendItem(mut, 'Float')
		per1 = tp.AppendItem(mut, 'Permutation')
		rec = tp.AppendItem(ga, 'Recombination')
		b2 = tp.AppendItem(rec, 'Boolean')
		i2 = tp.AppendItem(rec, 'Integer')
		d2 = tp.AppendItem(rec, 'Float')
		per2 = tp.AppendItem(rec, 'Permutation')
		psel = tp.AppendItem(ga, 'Parent Selection')
		ssel = tp.AppendItem(ga, 'Survivor Selection')

		list_ctrl = wx.ListCtrl(p2, 1, style=wx.LC_REPORT|wx.BORDER_SUNKEN)
#		list_ctrl = EditableListCtrl(self, size=(100, 100), style=wx.LC_REPORT|wx.BORDER_SUNKEN)
		list_ctrl.InsertColumn(0, 'Property')
		list_ctrl.InsertColumn(1, 'Value')
		list_ctrl.InsertStringItem(0, 'Population Size')
		list_ctrl.SetStringItem(0, 1, '100')
		list_ctrl.InsertStringItem(1, 'Maximum Iterations')
		list_ctrl.SetStringItem(1, 1, '50')

		plotgrid = wx.TextCtrl(p3, style=wx.TE_MULTILINE)

		stattx = wx.StaticText(p4, 1, label="Output goes here:")
		p4.SetBackgroundColour("White") # change to stattx.Set... later

		sz1.Add(tp, 1, wx.GROW)
		sz2.Add(list_ctrl, 1, wx.GROW)
		sz3.Add(plotgrid, 1, wx.GROW)
		sz4.Add(stattx, 1, wx.GROW)

		p1.SetSizer(sz1)
		p2.SetSizer(sz2)
		p3.SetSizer(sz3)
		p4.SetSizer(sz4)

		# initialize frame parameters
		self.SetSize((1024,768))
		self.SetTitle('Simple Toolbar')
		self.Centre()
		self.Show(True)

		self.InitUI()

	def InitUI(self):

		# initialize status bar
		self.CreateStatusBar()

		# define menu strings
		strNew = "Create a new document"
		strOpen = "Open a file"
		strSave = "Save the current file"
		strSaveAs = "Save the current file with a different name"
		strAppend = "Revert to a saved version of the file"
		strClose = "Close the current file"
		strExit = "Quit the program"
		strUndo = "Undo the last action"
		strRedo = "Redo the last undone action"
		strCut = "Cut the selection"
		strCopy = "Copy the selection"
		strPaste = "Paste the clipboard"
		strDelete = "Delete the selected text"
		strSelectAll = "Select the entire document"
		strProperties = "Configure the application"
		strHelp = "Open the help manual"
		strAbout = "About this application"

		# Setting up the menu.
		menuFile = wx.Menu()
		menuEdit = wx.Menu()
		menuView = wx.Menu()
		menuHelp = wx.Menu()

		# Populate menu
		menuFileNew = menuFile.Append(wx.ID_NEW, "&New", strNew)
		menuFileOpen = menuFile.Append(wx.ID_OPEN, "&Open...", strOpen)
		menuFile.AppendSeparator()
		menuFileSave = menuFile.Append(wx.ID_SAVE, "&Save", strSave)
		menuFileSaveAs = menuFile.Append(wx.ID_SAVEAS, "Save &As...", strSaveAs)
		menuFileRevert = menuFile.Append(wx.ID_REVERT, "&Revert", strAppend)
		menuFile.AppendSeparator()
		menuFileClose = menuFile.Append(wx.ID_CLOSE, "&Close", strClose)
		menuFileExit = menuFile.Append(wx.ID_EXIT, "&Quit", strExit)

		menuEditUndo = menuEdit.Append(wx.ID_UNDO, "&Undo", strUndo)
		menuEditRedo = menuEdit.Append(wx.ID_REDO, "&Redo", strRedo)
		menuEdit.AppendSeparator()
		menuEditCut = menuEdit.Append(wx.ID_CUT, "Cu&t", strCut)
		menuEditCopy = menuEdit.Append(wx.ID_COPY, "&Copy", strCopy)
		menuEditPaste = menuEdit.Append(wx.ID_PASTE, "&Paste", strPaste)
		menuEditDelete = menuEdit.Append(wx.ID_DELETE, "&Delete", strDelete)
		menuEdit.AppendSeparator()
		menuEditSelectAll = menuEdit.Append(wx.ID_SELECTALL, "Select &All", strSelectAll)
		menuEdit.AppendSeparator()
		menuEditProperties = menuEdit.Append(wx.ID_PROPERTIES, "Pr&eferences", strProperties)

		menuViewTree = menuView.Append(wx.ID_ANY, "Show tree", "Show tree", kind=wx.ITEM_CHECK)
		menuViewProperties = menuView.Append(wx.ID_ANY, "Show properties", "Show properties", kind=wx.ITEM_CHECK)
		menuViewScreen = menuView.Append(wx.ID_ANY, "Show screen", "Show screen", kind=wx.ITEM_CHECK)
		menuViewOutput = menuView.Append(wx.ID_ANY, "Show output", "Show output", kind=wx.ITEM_CHECK)

		menuHelpHelp = menuHelp.Append(wx.ID_HELP, "&Contents", strHelp)
		menuHelp.AppendSeparator()
		menuHelpAbout = menuHelp.Append(wx.ID_ABOUT, "&About", strAbout)

		# Creating the menubar.
		menuBar = wx.MenuBar()
		menuBar.Append(menuFile, "&File")
		menuBar.Append(menuEdit, "&Edit")
		menuBar.Append(menuView, "&View")
		menuBar.Append(menuHelp, "&Help")
		self.SetMenuBar(menuBar)

		# create toolbar
		toolbar = self.CreateToolBar()
		toolNew = toolbar.AddTool(wx.ID_ANY, wx.Bitmap('document-new.png'), shortHelpString='New', longHelpString=strNew)
		toolOpen = toolbar.AddTool(wx.ID_ANY, wx.Bitmap('document-open.png'), shortHelpString='Open', longHelpString=strOpen)
		toolSave = toolbar.AddTool(wx.ID_ANY, wx.Bitmap('document-save.png'), shortHelpString='Save', longHelpString=strSave)
		toolSaveAs = toolbar.AddTool(wx.ID_ANY, wx.Bitmap('document-save-as.png'), shortHelpString='Save As', longHelpString=strSaveAs)
		toolbar.AddSeparator()
		toolCopy = toolbar.AddTool(wx.ID_ANY, wx.Bitmap('edit-copy.png'), shortHelpString='Copy', longHelpString=strCopy)
		toolPaste = toolbar.AddTool(wx.ID_ANY, wx.Bitmap('edit-paste.png'), shortHelpString='Paste', longHelpString=strPaste)

		toolbar.Realize()

		# Set events.
#		self.Bind(wx.EVT_MENU, self.OnAbout, menuFileAbout)
#		self.Bind(wx.EVT_MENU, self.OnExit, menuFileExit)

class TreePanel(wx.Panel):
	def __init__(self, parent):
		wx.Panel.__init__(self, parent)
#		self.quote = wx.StaticText(self, label="Your quote:", pos=(20, 30))
		self.tree = wx.TreeCtrl(parent, 1, wx.DefaultPosition, (-1,-1), wx.TR_HIDE_ROOT|wx.TR_HAS_BUTTONS)
		root = self.tree.AddRoot('Genetic Algorithm')
		base = self.tree.AppendItem(root, 'Base Parameters')
		mut = self.tree.AppendItem(root, 'Mutation')
		rec = self.tree.AppendItem(root, 'Recombination')
		psel = self.tree.AppendItem(root, 'Parent Selection')
		ssel = self.tree.AppendItem(root, 'Survivor Selection')

def main():
	app = wx.App()
	MainWindow(None)
	app.MainLoop()

if __name__ == '__main__':
	main()


