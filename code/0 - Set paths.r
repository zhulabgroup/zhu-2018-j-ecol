# set paths
file.symlink("/data/ZHULAB/FIA/Mycorrhiza/Input/", "Data")
file.symlink("/data/ZHULAB/FIA/Mycorrhiza/Models/", "Models")

if (!file.exists("Figures")) {
  dir.create("Figures")
}
