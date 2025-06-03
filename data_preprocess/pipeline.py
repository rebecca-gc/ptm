import download_glycosylation
import download_s_nitrosylation
import download_acetylation
import download_methylation
import merge
import negatives
import class_generator


def main():
    download_glycosylation.main()
    download_s_nitrosylation.main()
    download_acetylation.main()
    download_methylation.main()
    merge.main()
    negatives.main()
    class_generator.main()

    print("\nEverything worked! :)\n")


if __name__ == '__main__':
    main()
