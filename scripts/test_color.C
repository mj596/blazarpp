void test_color(Int_t r, Int_t g, Int_t b)
{
  TApplication::NeedGraphicsLibs();
  gApplication->InitializeGraphics();
  TString color = TString::Format("#%06x", TColor::RGB2Pixel(r, g, b));
  cout << color << endl;
  TColor col;
  col.SetRGB((Float_t)r/255.0, (Float_t)g/255.0, (Float_t)b/255.0);
  cout << col.AsHexString() << endl;
}
