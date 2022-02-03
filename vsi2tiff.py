import argparse
import bioformats
import javabridge
import os
import re
import subprocess
from utils import run_command


def _init_logger():
    """This is so that Javabridge doesn't spill out a lot of DEBUG messages
    during runtime.
    From CellProfiler/python-bioformats.
    """
    rootLoggerName = javabridge.get_static_field("org/slf4j/Logger",
                                         "ROOT_LOGGER_NAME",
                                         "Ljava/lang/String;")

    rootLogger = javabridge.static_call("org/slf4j/LoggerFactory",
                                "getLogger",
                                "(Ljava/lang/String;)Lorg/slf4j/Logger;",
                                rootLoggerName)

    logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",
                                   "WARN",
                                   "Lch/qos/logback/classic/Level;")

    javabridge.call(rootLogger,
            "setLevel",
            "(Lch/qos/logback/classic/Level;)V",
            logLevel)

parser = argparse.ArgumentParser()

parser.add_argument('--input_vsi', type=str, required=True)
#parser.add_argument('--sample_name', type=str, required=True)
parser.add_argument('--output_dir', type=str, required=True)
parser.add_argument('--bftools_dir', type=str, default='/mnt/data2/shared/Lymphomoid-IF-software/bftools')

args = parser.parse_args()

if not os.path.splitext(args.input_vsi):
    raise ValueError('The path from --input_vsi should be the path to the .vsi file.')

parent_path, vsi_name = os.path.split(args.input_vsi)
vsi_name = os.path.splitext(vsi_name)[0]
sample_name = os.path.split(parent_path)[1]

javabridge.start_vm(class_path=bioformats.JARS)
try:
    _init_logger()

    omexml = bioformats.get_omexml_metadata(args.input_vsi)
    o = bioformats.OMEXML(omexml)

    regex = "[0-9]{2}x_[0-9]{2}"

    count = 1
    for i in range(o.image_count):
        if re.match(regex, o.image(i).Name):
            output_path = os.path.join(args.output_dir, 'Tiff', f'{sample_name}_{vsi_name}_acq{count:02d}.ome.tif')
            if os.path.exists(output_path):
                os.remove(output_path)
            run_command(f"{os.path.join(args.bftools_dir, 'bfconvert')} -option BF_MAXMEM 128g -series 12  -bigtiff -pyramid-resolutions 1 {args.input_vsi} {output_path}")
            print(f'Acquisition {count} exported to {output_path}')
            print()
            count += 1
except Exception as e:
    javabridge.kill_vm()
    raise e
javabridge.kill_vm()