### README

#### Notice
This project is provided *as-is*, with no warranties whatsoever, no guarantees, etc. This constitutes a class project, and depth for project scope and

#### Project Scope

The aim of this project is to balance an inherently unstable system with motion adjustments provided by an external source. Specifically, a DC motor is used to provide motion adjustments to balance an inverted pendulum.


#### Project Layout
```
.
├── cad
├── figs
│   ├── drawing_tree.drawio
│   ├── drawing_tree.drawio.png
│   ├── example_simulation.png
│   └── qr_code_to_repo.png
├── poster
│   ├── poster_here.pdf
│   ├── poster.kra
│   ├── poster.kra~
│   ├── poster.pdf
│   ├── poster.png
│   ├── poster.png.kra
│   └── poster.tex
├── README.md
├── schem
├── sim
│   ├── requirements.txt
│   └── simulation.py
├── sketch
│   └── sketch.ino
└── summary.txt

7 directories, 16 files
```

#### Quickstart

```bash 
git clone https://github.com/A-Person7/MENG-451-final-project
cd sim
pip3 install -r requirements.txt # If you need to, may want to setup a venv first if your system is externally managed
python3 simulation.py
```

#### Poster
To compile the poster to a PDF/similar with `latexmk` installed, 
```bash
cd poster
latexmk -pvc -xelatex poster.tex
```
