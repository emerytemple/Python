#!/usr/bin/python

import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import wx

class Example(wx.Frame):

	def __init__(self, *args, **kwargs):
		super(Example, self).__init__(*args, **kwargs)

		self.InitUI()

	def InitUI(self):

		pnl = wx.Panel(self)

		hbox_main = wx.BoxSizer(wx.HORIZONTAL)
		vbox = wx.BoxSizer(wx.VERTICAL)
		hbox1 = wx.BoxSizer(wx.HORIZONTAL)
		hbox2 = wx.BoxSizer(wx.HORIZONTAL)
		hbox3 = wx.BoxSizer(wx.HORIZONTAL)
		hbox4 = wx.BoxSizer(wx.HORIZONTAL)
		hbox5 = wx.BoxSizer(wx.HORIZONTAL)
		hbox6 = wx.BoxSizer(wx.HORIZONTAL)
		hbox7 = wx.BoxSizer(wx.HORIZONTAL)
		hbox8 = wx.BoxSizer(wx.HORIZONTAL)
		hbox9 = wx.BoxSizer(wx.HORIZONTAL)
		hbox10 = wx.BoxSizer(wx.HORIZONTAL)
		hbox11 = wx.BoxSizer(wx.HORIZONTAL)
		hbox12 = wx.BoxSizer(wx.HORIZONTAL)

		st1 = wx.StaticText(pnl, label='n1')
		st2 = wx.StaticText(pnl, label='n2')
		st3 = wx.StaticText(pnl, label='peipc')
		st4 = wx.StaticText(pnl, label='xp')
		st5 = wx.StaticText(pnl, label='rrre')
		st6 = wx.StaticText(pnl, label='pht')
		st7 = wx.StaticText(pnl, label='r')
		st8 = wx.StaticText(pnl, label='te')
		st9 = wx.StaticText(pnl, label='papc')
		st10 = wx.StaticText(pnl, label='g')
		st11 = wx.StaticText(pnl, label='gamma')
		st12 = wx.StaticText(pnl, label='trunc')

		self.tc1 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc2 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc3 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc4 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc5 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc6 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc7 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc8 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc9 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc10 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc11 = wx.TextCtrl(pnl, size=(180, -1))
		self.tc12 = wx.TextCtrl(pnl, size=(180, -1))

		#self.tc = wx.TextCtrl(pnl, style=wx.TE_MULTILINE)
		button_send = wx.Button(pnl, label='Execute')

		self.statusbar = self.CreateStatusBar()
		self.statusbar.SetStatusText("Hola")

		hbox1.Add(st1, flag=wx.LEFT, border=10)
		hbox1.Add(self.tc1, flag=wx.LEFT, border=10)
		hbox2.Add(st2, flag=wx.LEFT, border=10)
		hbox2.Add(self.tc2, flag=wx.LEFT, border=10)
		hbox3.Add(st3, flag=wx.LEFT, border=10)
		hbox3.Add(self.tc3, flag=wx.LEFT, border=10)
		hbox4.Add(st4, flag=wx.LEFT, border=10)
		hbox4.Add(self.tc4, flag=wx.LEFT, border=10)
		hbox5.Add(st5, flag=wx.LEFT, border=10)
		hbox5.Add(self.tc5, flag=wx.LEFT, border=10)
		hbox6.Add(st6, flag=wx.LEFT, border=10)
		hbox6.Add(self.tc6, flag=wx.LEFT, border=10)
		hbox7.Add(st7, flag=wx.LEFT, border=10)
		hbox7.Add(self.tc7, flag=wx.LEFT, border=10)
		hbox8.Add(st8, flag=wx.LEFT, border=10)
		hbox8.Add(self.tc8, flag=wx.LEFT, border=10)
		hbox9.Add(st9, flag=wx.LEFT, border=10)
		hbox9.Add(self.tc9, flag=wx.LEFT, border=10)
		hbox10.Add(st10, flag=wx.LEFT, border=10)
		hbox10.Add(self.tc10, flag=wx.LEFT, border=10)
		hbox11.Add(st11, flag=wx.LEFT, border=10)
		hbox11.Add(self.tc11, flag=wx.LEFT, border=10)
		hbox12.Add(st12, flag=wx.LEFT, border=10)
		hbox12.Add(self.tc12, flag=wx.LEFT, border=10)

		vbox.Add(hbox1, flag=wx.TOP, border=10)
		vbox.Add(hbox2, flag=wx.TOP, border=10)
		vbox.Add(hbox3, flag=wx.TOP, border=10)
		vbox.Add(hbox4, flag=wx.TOP, border=10)
		vbox.Add(hbox5, flag=wx.TOP, border=10)
		vbox.Add(hbox6, flag=wx.TOP, border=10)
		vbox.Add(hbox7, flag=wx.TOP, border=10)
		vbox.Add(hbox8, flag=wx.TOP, border=10)
		vbox.Add(hbox9, flag=wx.TOP, border=10)
		vbox.Add(hbox10, flag=wx.TOP, border=10)
		vbox.Add(hbox11, flag=wx.TOP, border=10)
		vbox.Add(hbox12, flag=wx.TOP, border=10)
		vbox.Add(button_send, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=20)

		self.figure = matplotlib.figure.Figure()
		self.axes = self.figure.add_subplot(111)
		t = [0,1,2,3,4,5,6,7,8,9]
		s = [0,1,0,1,0,2,1,2,1,0]
		self.y_max = 10
		self.axes.plot(t,s)
		self.canvas = FigureCanvas(pnl, -1, self.figure)

		hbox_main.Add(self.canvas, proportion=1, flag=wx.EXPAND | wx.TOP | wx.RIGHT | wx.LEFT, border=15)
		hbox_main.Add(vbox, flag=wx.TOP, border=10)

		# Callbacks
		self.Bind(wx.EVT_BUTTON, self.OnSend, button_send)

		pnl.SetSizer(hbox_main)
		self.SetSize((400, 420))
		self.SetTitle('Aerospike Programme')
		self.Centre()
		self.Show()

	def OnSend(self,event):

		aerospike = Curve()
		xx, yy = aerospike.calc_curve()
		print xx, yy

class Curve:

	def __init__(self):
		self.n1 = 50
		self.n2 = 500
		self.peipc = 0.0994474
		self.xp = 12.0
		self.rrre = 0.3
		self.pht = 1.570796
		self.r = 1718.0
		self.te = 3000.0
		self.papc = 0.01
		self.g = 32.17405
		self.gamma = 1.2
		self.trunc = 1.0

		self.tol = 1e-14
		self.gp1 = self.gamma + 1
		self.gm1 = self.gamma - 1
		self.exp = self.gp1/(2*self.gm1)
		self.c1 = 2/self.gp1
		self.c2 = self.gm1/self.gp1
		self.c3 = self.gp1/self.gm1
		self.c4 = (1-self.gamma)/self.gamma
		self.c5 = 2/self.gm1
		self.c6 = self.gamma/self.gm1

	def mach_func(self, M):
		fme = ((M*self.xp)-((self.c1+(self.c2*M*M))**self.exp))
		return fme

	def find_mach(self):
		# bisection
		for i in range(2, 51): 
			x = self.mach_func(i)

			if x < 0.0:
				a = i-1
				b = i
				break

		u = self.mach_func(a)
		v = self.mach_func(b)

		for n in range(self.n1):
			c = (a + b)/2.0
			w = self.mach_func(c)
			if math.fabs(w*u) < self.tol:
				ans = c
				break
			elif w*u < 0:
				b = c
				v = w
			else:
				a = c
				u = w
		return ans

	def find_mu(self, M):
		fmu = ((math.sqrt(self.c3)*math.atan(math.sqrt(self.c2*((M*M)-1.0))))-math.atan(math.sqrt((M*M)-1.0)))
		return fmu

	def calc_curve(self):
		rm = self.find_mach()
		ve = self.find_mu(rm)
		delta = (math.pi/2.0)-ve
		sd = math.sin(delta)
		ht = (self.xp-math.sqrt(self.xp*(self.xp-sd)))/(self.xp*sd)
		cfo = self.gamma*rm*(self.c1**self.exp)/math.sqrt(1.0+0.5*self.gm1*rm*rm)
		pcpt = (0.5*self.gp1)**self.c6
		vt = math.sqrt((self.gamma*self.r*self.te*(1.0+0.5*self.gm1*rm*rm)/(0.5*self.gp1))*self.g)
		spim = 1.0+((1.0-pcpt*self.papc)/self.gamma)
		
		drm = (rm - 1.0)/self.n1

		xm = []
		rxre = []
		xxre = []
		pxpc = []
		sumcg = []
		sumim = []
		sumva = []

		xm.append(1.0)
		rxre.append(1.0-ht*sd)
		xxre.append((-ht)*math.cos(delta))
		pxpc.append((1.0+0.5*self.gm1*xm[0]*xm[0])**(-self.c6))
		sumcg.append((self.c1**self.c6)*self.gp1*sd)
		sumim.append(vt*sd*spim/self.g)
		sumva.append(vt*sd*(1.0+1.0/self.gamma)/self.g)

#		print xm[0], rxre[0], xxre[0], pxpc[0], sumcg[0], sumim[0], sumva[0]

		for k in range(1,self.n1+1):
			xm.append(xm[k-1] + drm)
			pro = pxpc[k-1]
			rxo = rxre[k-1]
			vx = self.find_mu(xm[k])
			y = 1.0/xm[k]
			ux = math.atan(y/math.sqrt(1.0-y*y))
			xm2 = xm[k]*xm[k]
			rxre.append(math.sqrt(math.fabs(1.0-((self.c1*(1.0+0.5*self.gm1*xm2))**self.exp)*math.sin(ve-vx+ux)/self.xp)))
			xxre.append((1.0-rxre[k])*math.cos(ve-vx+ux)/math.sin(ve-vx+ux))
			pxpc.append((1.0+0.5*self.gm1*xm2)**(-self.c6))

			rr = rxo*rxo-rxre[k]*rxre[k]

			sumcg.append(sumcg[k-1]+0.5*self.xp*(pro+pxpc[k])*rr)
			co = pcpt*vt*self.xp/(self.g*self.gamma)
			sumim.append(sumim[k-1]+0.5*co*(pro+pxpc[k]-2.0*self.papc)*rr)
			sumva.append(sumva[k-1]+0.5*co*(pro+pxpc[k])*rr)
#			print xm[k], rxre[k], xxre[k], pxpc[k], sumcg[k], sumim[k], sumva[k]
		return xm, rxre

	def truncate(self, x, y):
		xx = []
		yy = []
		xlen = x[self.n1] - x[0]
		tlen = self.trunc*xlen
		tmark = x[0] + tlen
		for i in range(self.n1+1):
			if x[i] < tmark:
				xx.append(x[i])
				yy.append(y[i])
		return xx, yy

def main():

	ex = wx.App()
	Example(None)
	ex.MainLoop()

if __name__ == '__main__':
	main()
