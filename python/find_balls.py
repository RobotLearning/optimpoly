import re
import tensorflow as tf
import os
import scipy
import numpy as np

model_name = '/home/okoc/Dropbox/capture_train/trained/lin_lr'


def find_balls(img_path, ranges=None, cams='all'):
    '''
    Find the ball pixels for every image in the image path
    that fall in the range given for given cameras

    Returns the pixels as a Nx2 np array
    '''
    fexp = re.compile(r'([\w]+)\.(jpg)')
    examples_list = filter(fexp.search, os.listdir(img_path))
    imgs = []

    if cams is not 'all':
        s = str(cams[0])
        for i in range(len(cams)-1):
            s = s + r'|' + str(cams[i+1])
        fexp2 = re.compile('c(' + s + ')*')
        examples_list = filter(fexp2.search, examples_list)

    ex_list = []
    if ranges is not None:
        for idx, ex in enumerate(examples_list):
            s = ex.split('.')
            ss = s[0].split('_')
            num = int(ss[1])
            if num >= ranges[0] and num <= ranges[1]:
                ex_list.append(ex)
    else:
        ex_list = examples_list

    with tf.Session() as sess:
        new_saver = tf.train.import_meta_graph(
            '{}.meta'.format(model_name))
        new_saver.restore(sess, '{}'.format(model_name))
        conv_params = sess.run('conv2d/kernel:0')
        bias_params = sess.run('conv2d/bias:0')
        print("log reg params: {}".format(conv_params[0, 0, :, 1]))
        print("bias: {}".format(bias_params))

    pixels_ball = np.zeros((len(ex_list), 2))
    for i, example in enumerate(ex_list):
        img = scipy.ndimage.imread(os.path.join(img_path, example))
        tf.reset_default_graph()

        prob = tf.get_collection("probabilities")[0]
        graph_input = img/255.0
        prob_img = prob.eval(feed_dict={x: graph_input})
        idx_max = prob_img.argmax()
        idx_row = idx_max / prob_img.shape[0]
        idx_col = idx_max - (idx_row * prob_img.shape[0])
        pixels_ball[i, 0] = idx_row
        pixels_ball[i, 1] = idx_col
    return pixels_ball
