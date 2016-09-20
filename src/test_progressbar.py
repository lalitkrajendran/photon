import time
from progressbar import ProgressBar

pbar = ProgressBar(maxval=5).start()

for i in range(0,6):
    pbar.update(i)
    time.sleep(1)
pbar.finish()