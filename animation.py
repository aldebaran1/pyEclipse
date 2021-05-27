import cv2
import os
from argparse import ArgumentParser
import numpy as np


def lapse(folder: str = None,
          fps: float = None,
          name: str = None,
          gamma: float = None,
          scale: float = None,
          extension: str = '.png'
          ):
    '''convert sequential png files to mp4'''

    invgamma = 1 / gamma
    table = np.array([((i / 255.0) ** invgamma) * 255 for i in np.arange(0, 256)]).astype('uint8')

    fourcc = 0x7634706d
    images = sorted([img for img in os.listdir(folder) if img.endswith(extension)])
    frame = cv2.imread(os.path.join(folder, images[0]))

#    if height is None or width is None:
    
    if name is None:
        name = folder + 'movie.mp4'
#        name = os.path.split(folder)[1]+'.mp4'
    height, width, layers = frame.shape
    if scale is not None:
        width = int(width * scale)
        height = int(height * scale)
    video = cv2.VideoWriter(name, fourcc, fps, (width, height))

    for image in images:
        im = cv2.imread(os.path.join(folder, image))
        res = cv2.resize(cv2.LUT(im, table), (width, height), interpolation=cv2.INTER_AREA)
        video.write(res)

    cv2.destroyAllWindows()
    video.release()


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('folder', type=str, help='folder containing images')
    p.add_argument('-s', '--scale', type=float, help='rescale the image - relative factor', default=None)
    p.add_argument('-f', '--fps', type=float, help='frames per second', default=5)
    p.add_argument('-n', '--name', type=str, help='output name with .mp4 ext', default=None)
    p.add_argument('-g', '--gamma', type=float, default=1.0)
    p.add_argument('-e', '--sufix', type=str, default='.png', help="File extenion. Default: .png")

    P = p.parse_args()

    lapse(folder=P.folder, scale=P.scale, fps=P.fps, name=P.name, gamma=P.gamma, extension=P.sufix)
