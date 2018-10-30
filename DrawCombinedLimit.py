import ROOT
import sys


CMSfile = ROOT.TFile("LimitcontoursAxial_observed_hist_bothCMS.root","OPEN")
combinedfile = ROOT.TFile("LimitcontoursAxial_observed_hist_bothcombined.root","OPEN")

cms = ROOT.TLatex(0.12, 0.89, '#splitline{CMS #scale[0.8]{#it{simulation}}}{#scale[0.8]{work in progress}}'  )
cms.SetNDC()
cms.SetTextSize(0.04)
ROOT.gStyle.SetPalette(1)
lumi = ROOT.TLatex(0.7, 0.91, '35.9 fb^{-1} (13 TeV)'  )
lumi.SetNDC()
lumi.SetTextSize(0.03)

ROOT.gStyle.SetPadRightMargin(0.15)
w=1200
h=800
canvas = ROOT.TCanvas()
canvas.cd()
canvas.SetCanvasSize(w,h)

expGraph = combinedfile.Get("MCgraph_expected")
ROOT.gPad.SetLogz()
expGraph.Draw("colz")
expGraph.GetHistogram().GetXaxis().SetTitle("m_{med} [GeV/c^{2}]")
expGraph.GetHistogram().GetXaxis().SetTitleSize(0.05)
expGraph.GetHistogram().GetXaxis().SetTitleOffset(0.8)

expGraph.GetHistogram().GetYaxis().SetTitle("m_{DM} [GeV/c^{2}]")
expGraph.GetHistogram().GetYaxis().SetTitleSize(0.05)
expGraph.GetHistogram().GetYaxis().SetTitleOffset(0.8)

expGraph.GetHistogram().SetMinimum(0.01)
expGraph.GetHistogram().GetZaxis().SetTitle("#sigma_{95% CL}^{expected}/#sigma_{th}")
expGraph.GetHistogram().GetZaxis().SetTitleSize(0.05)
expGraph.GetHistogram().GetZaxis().SetTitleOffset(0.8)

expGraph.SetTitle("Exclusion limits at 95% CL")

expContcombined = combinedfile.Get("expected")
expContcombined.Draw("same")
expContcombined.SetLineColor(ROOT.kRed)

expContUp1 = combinedfile.Get("up1")
expContUp1.Draw("same")
expContUp1.SetLineColor(ROOT.kRed)
expContUp1.SetLineStyle(5)


expContDown1 = combinedfile.Get("down1")
expContDown1.Draw("same")
expContDown1.SetLineColor(ROOT.kRed)
expContDown1.SetLineStyle(5)

expContCMS = CMSfile.Get("expected")
expContCMS.Draw("same")

expContCMSUp1 = CMSfile.Get("up1")
expContCMSUp1.Draw("same")
expContCMSUp1.SetLineStyle(5)


expContCMSDown1 = CMSfile.Get("down1")
expContCMSDown1.Draw("same")
expContCMSDown1.SetLineStyle(5)

points=CMSfile.Get("points")
points.Draw("*")

legend = ROOT.TLegend(0.1,0.7,0.48,0.9);
legend.AddEntry(expContcombined,"Combined median expected","l");
legend.AddEntry(expContUp1,"Combined median expected #pm 1 #sigma","l");
legend.AddEntry(expContCMS,"Median expected","l");
legend.AddEntry(expContCMSUp1,"Median expected #pm 1 #sigma","l");
legend.Draw();

canvas.SaveAs("combinedLimitAxial.pdf")