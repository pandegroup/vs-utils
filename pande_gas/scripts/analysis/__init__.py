import pandas as pd


def get_scores(filename):
    """
    Get scores and dataset divisions.
    """
    df = pd.read_table(filename)
    scores = {}
    ref_idx = 1
    new_idx = 27
    datasets = {'PCBA': [], 'MUV': [], 'TOX': [], 'DUDE': []}
    print df.values[ref_idx][0], 'VS.', df.values[new_idx][0]  # print scores
    for name, ref_score, new_score in zip(
            df.columns[1:], df.values[ref_idx][1:], df.values[new_idx][1:]):
        score = new_score - ref_score
        if name.startswith('PCBA'):
            name = name.split('PCBA-AID')[-1]
            datasets['PCBA'].append(name)
        elif name.startswith('MUV'):
            name = name.split('MUV-')[-1]
            datasets['MUV'].append(name)
        elif name.startswith('TOX'):
            name = name.split('-')
            name.pop()
            name.pop(0)
            name = '_'.join(name)
            datasets['TOX'].append(name)
        elif name.startswith('DUDE'):
            name = name.split('DUDE-')[-1]
            datasets['DUDE'].append(name)
        else:
            raise ValueError(name)
        scores[name] = score

    # sanity checks
    total = 0
    for key in datasets:
        total += len(datasets[key])
    assert total == 259, total

    return scores, datasets
