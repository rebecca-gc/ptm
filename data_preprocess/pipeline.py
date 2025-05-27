import download_glycosylation
import download_s_nitrosylation
import download_acetylation
import merge
import negatives


def main():
    download_glycosylation.main()
    download_s_nitrosylation.main()
    download_acetylation.main()
    merge.main()
    negatives.main()

    print("\nEverything worked! :)\n")


if __name__ == '__main__':
    main()
