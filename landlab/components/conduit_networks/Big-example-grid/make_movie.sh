#!/bin/sh

mencoder "mf://*.png" -mf type=png:fps=18 -ovc lavc -o animation.avi
