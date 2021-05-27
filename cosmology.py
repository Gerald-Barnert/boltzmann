import camb
from classy import Class
import numpy as np
import matplotlib.pyplot as plt
import final

camb_code = final.run("camb", 70, 0.04, 0.24)
class_code = final.run("class", 70, 0.04, 0.24)