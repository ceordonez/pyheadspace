from scr.process_data import process_data
from scr.read_data import read_config, read_data
from scr.write_data import write_data


def main():

    cfg = read_config("config.yml")
    data = read_data(cfg)
    processdata = process_data(cfg, data)
    write_data(cfg, processdata)


if __name__ == "__main__":
    main()
