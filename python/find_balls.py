import re
import tensorflow as tf
import os
import scipy
import numpy as np
from matplotlib import pyplot as plt

#model_name = os.environ['HOME'] + '/vision/Vision_Experiments/sebastian/trained/lin_lr'
model_name = os.environ['HOME'] + \
    '/Dropbox/capture_train/trained/lin_lr'  # new version


def find_balls(img_path, ranges, cams='all', debug=False, prefix='c'):
    '''
    Find the ball pixels for every image in the image path
    that fall in the range given for given cameras (e.g. 100-200 for cams 0 and 1)
    Returns the pixels as a dictionary from image number to pixels (at least 4D)
    '''
    fexp = re.compile(r'([\w]+)\.(jpg)')
    examples_list = filter(fexp.search, os.listdir(img_path))
    imgs = []
    #print examples_list

    if cams is not 'all':
        s = str(cams[0])
        for i in range(len(cams)-1):
            s = s + r'|' + str(cams[i+1])
        fexp2 = re.compile(prefix + '([' + s + '])')
        examples_list = filter(fexp2.search, examples_list)

    ex_dict = dict()
    for idx, ex in enumerate(examples_list):
        s = ex.split('.')
        ss = s[0].split('_')
        num = int(ss[1])
        if num >= ranges[0] and num <= ranges[1]:
            try:
                ex_dict[num].append(ex)
                ex_dict[num].sort()
            except KeyError:
                ex_dict[num] = [ex]

    # print(ex_dict)
    tf.reset_default_graph()
    with tf.Session() as sess:
        new_saver = tf.train.import_meta_graph(
            '{}.meta'.format(model_name))
        new_saver.restore(sess, '{}'.format(model_name))
        conv_params = sess.run('conv2d/kernel:0')
        bias_params = sess.run('conv2d/bias:0')
        print("log reg params: {}".format(conv_params[0, 0, :, 1]))
        print("bias: {}".format(bias_params))

        pixels_ball = dict()
        for i, examples in ex_dict.iteritems():
            if len(examples) >= 2:
                idxs = []
                prob_imgs = []
                imgs = []
                for j, example in enumerate(examples):
                    imgs.append(scipy.ndimage.imread(
                        os.path.join(img_path, example)))
                    prob = tf.get_collection("probabilities")[0]
                    x = tf.get_collection("input")[0]
                    graph_input = imgs[-1]/255.0
                    graph_input = graph_input[np.newaxis]  # one by one
                    prob_img_tf = prob.eval(feed_dict={x: graph_input})
                    # print(prob_img.shape)
                    prob_imgs.append(np.squeeze(prob_img_tf[0, :, :, 1]))
                    idx_max = prob_imgs[-1].argmax()
                    idx_row = idx_max / prob_imgs[-1].shape[1]  # row
                    idx_col = idx_max - \
                        (idx_row * prob_imgs[-1].shape[1])  # col
                    idxs.append(idx_col)
                    idxs.append(idx_row)
                pixels_ball[i] = np.array(idxs)
                if debug:
                    # print(prob_img.shape)
                    print(pixels_ball[i])
                    fig = plt.figure(figsize=(32, 12))
                    for j, prob_img in enumerate(prob_imgs):
                        a = fig.add_subplot(1, 2*len(prob_imgs), 2*j+1)
                        a.imshow(imgs[j])
                        b = fig.add_subplot(1, 2*len(prob_imgs), 2*j+2)
                        b.imshow(prob_imgs[j], cmap='gray')
                    plt.show()
                    raw_input("Press Enter to continue...")

            if i % 100 == 0:
                print('Evaluated prob for image', i)

    return pixels_ball
