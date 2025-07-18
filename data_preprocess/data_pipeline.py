import os
import download_all
import merge
import venn_diagramms
import negatives
import class_generator


def main():
    #download_all.main()
    ptms_dir = 'data/ptms'
    for dir in os.listdir(ptms_dir):
        if os.path.isdir(dir):
            dir_path = os.path.join(ptms_dir, dir)
            databases_path = os.path.join(dir_path, 'databases')
            merge.main(databases_path, dir_path)
    venn_diagramms.main(ptms_dir)
    negatives.main()
    for dir in os.listdir(ptms_dir):
        dir_path = os.path.join(ptms_dir, dir)
        if os.path.isdir(dir_path):
            class_generator.generator(os.path.join(dir_path, 'merged3000.fasta'),f'data/no_ptm/filtered3000_no_{dir}.fasta',dir_path)
            databases_path = os.path.join(dir_path, 'databases')
            for db in os.listdir(databases_path):
                class_generator.generator(os.path.join(databases_path, db),f'data/no_ptm/filtered3000_no_{dir}.fasta',os.path.join(dir_path, 'db_seqs_classes'),db=f'{dir}-{db.split(".")[0]}_')


    print('\nEverything worked! :)\n')


if __name__ == '__main__':
    main()
