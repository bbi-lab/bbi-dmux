#
# run_info.py reads an Illumina [rR]unParameters.xml file
#             from in an Illumina run directory.
#
# Notes:
#   o  the master file is bbi-dmux-data/illumina_run_parameters/src/run_info.py
#      Make changes to the master file and distribute
#      to all other pipelines that use it in order to
#      maintain consistency.
#   o  this file is used in the following pipelines
#         o  bbi-dmux: sci-RNA-seq demultiplexing pipeline
#         o  sciatac_pipeline: sci-ATAC-seq processing pipeline
#         o  bbi-sciatac-demux: sci-ATAC-seq demultiplexing pipeline
#   o  bbi-dmux/data/illumina_run_parameters has a collection of
#      [rR]unParameters.xml files for testing run_info.py. The
#      script bbi-dmux/data/illumina_run_parameters/run_test.sh
#      runs read_run_info.py on the example runParameters.xml files.
#

import os
import xml.etree.ElementTree as ET
import glob
import gzip
import re
from collections import OrderedDict
from typing import Union
from xml.etree.ElementTree import ElementTree, Element

#
# **  application_name identifiers **
# The <ApplicationName> element in RunParameters.xml
# is compared to a regular expression consisting of
# an RE alternation of the application_name_dict keys.
# The corresponding dictionary value is assigned to
# instrument_name in get_application_info().
#
# Notes:
#   o  the application_name string identifies the
#      instrument_type, which determines the value of
#      of reverse_complement_i5. reverse_complement_i5
#      is used in barcode_correct_sciatac.py to use
#      the correct reverse_complement_i5 orientation.
#      That is, the application_name and
#      SEQUENCERS_P5_RC_MAP are required for
#      successful processing.
#   o  the RE strings must distinguish the instruments
#   o  the <ApplicationVersion> value qualifies
#      the HiSeq instrument_name to HiSeq3000 or
#      HiSeq4000 in get_application_info().
#   o  number of lanes
#        MiSeq - 1 lane cells
#        HiSeg machines had 8 lanes in High Output Mode and 2 lane cells in Rapid Mode (only for HiSeq 2000/2500 or other too??)
#        NextSeq - 4 lanes
#        NextSeq 1000/2000 - 1 lane
#        X Ten - 8 lane cells
#        NovaSeq S1, S2 - 2 lanes flowcells
#        NovaSeq S3, S4 - 4 lanes flowcells
#        NovaSeqX use SampleTube ConsumableInfo to determine number of lanes (I don't know whether this is definitive.)
#      note: the NextSeq and NovaSeq have 'fluidically' one lane; except for
#            NovaSeq XP protocol, which allows for distinct samples on lanes
#   o  few run_status values are required by the
#      pipelines. These include
#        o  reverse_complement_i5
#        o  p5_index_length
#        o  p7_index_length
#        o  lanes  (used in Andrew's sciatac_pipeline code)
#      Additional values may be recorded in the log
#      files.


version = '20240530.1'


application_name_dict = {
    'NextSeq Control': 'NextSeq500',
    'NextSeq 1000/2000': 'NextSeq2000',
    'MiSeq': 'MiSeq',
    'NovaSeqXSeries': 'NovaSeqX',
    'NovaSeq': 'NovaSeq',
    'HiSeq': 'HiSeq' }


def open_file(f, mode='rt'):
    if f.endswith('.gz'):
        return gzip.open(f, mode)
    else:
        return open(f, mode)

#
# Convert etree to dictionary.
# This is unused at this time but it may be useful
# in the future.
#
def etree_to_dict(root: Union[ElementTree, Element], include_root_tag=False):
    root = root.getroot() if isinstance(root, ElementTree) else root
    result = OrderedDict()
    if len(root) > 1 and len({child.tag for child in root}) == 1:
        result[next(iter(root)).tag] = [etree_to_dict(child) for child in root]
    else:
        for child in root:
            result[child.tag] = etree_to_dict(child) if len(list(child)) > 0 else (child.text or "")
    result.update(('@' + k, v) for k, v in root.attrib.items())
    return {root.tag: result} if include_root_tag else result


def get_application_info( tree ):
    """
    Identify the sequencing machine model using the RunParameters.xml file contents.
    Args:
      xml tree
    Returns:
      instrument_model: string with instrument name (known names are in re_models below)'
      application_version: string with run application version
    """

    # Search for an <ApplicationName> tag.
    application_name = None
    # most machines store the machine name string in the tag 'ApplicationName'
    for application_name in tree.getroot().iter( 'ApplicationName' ):
        application_name = application_name.text
        break

    # NovaSeq stores the machine name string in the tag 'Application'
    # so search for <Application> if application_name is None.
    if( application_name == None ):
        for application_name in tree.getroot().iter( 'Application' ):
            application_name = application_name.text
            break

    if( application_name == None ):
        raise ValueError( 'Unable to find Application* element in BCL RunParameters.xml' )

    application_version = None
    for application_version in tree.getroot().iter( 'ApplicationVersion' ):
        application_version = application_version.text
        break

    # NovaSeqX stores the application version in the tag SystemSuiteVersion
    if( application_version == None ):
       for application_version in tree.getroot().iter( 'SystemSuiteVersion' ):
           application_version = application_version.text
           break

    if( application_version == None ):
        raise ValueError( 'ApplicationVersion element missing in BCL RunParameters.xml' )

    re_models = '|'.join( application_name_dict.keys() )
    re_pattern = '(%s)' % re_models
    mobj = re.match( re_pattern, application_name )
    if( mobj == None ):
        raise ValueError( 'unrecognized ApplicationName in RunParameters.xml file' )
    instrument_model = application_name_dict[mobj.group( 1 )]

    # Distinguish between HiSeq models 3000 and 4000 using Andrew's(?) method.
    # Note: the p5 index orientations differ between these two models.
    if( instrument_model == 'HiSeq' ):
        application_major_version = int(application_version.split('.')[0])
        if application_major_version > 2:
            instrument_model = 'HiSeq4000'
        else:
            instrument_model = 'HiSeq3000'

    return( instrument_model, application_version )


def get_run_info( flow_cell_path, pipeline_type='RNA-seq' ):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        flow_cell_path (str): Path to BCL directory for run.
        pipeline_type (str): Processing pipeline type. Either 'RNA-seq' or 'ATAC-seq'.
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """

    bcl_run_info = os.path.join( flow_cell_path, '[rR]unParameters.xml*' )
    bcl_run_info = glob.glob( bcl_run_info )
    if( not bcl_run_info ):
        raise ValueError('BCL RunParameters.xml not found for specified flowcell: %s' % bcl_run_info)
    else:
        bcl_run_info = bcl_run_info[0]

    # Set up a few nodes for parsing
    tree = ET.parse(open_file(bcl_run_info))

    ( instrument_model, application_version ) = get_application_info( tree )

    if( instrument_model == 'NextSeq500' ):
        run_stats = get_run_info_nextseq500( instrument_model, application_version, tree )
    elif( instrument_model == 'NextSeq2000' ):
        run_stats = get_run_info_nextseq2000( instrument_model, application_version, tree )
    elif( instrument_model == 'MiSeq' ):
        run_stats = get_run_info_miseq( instrument_model, application_version, tree )
    elif( instrument_model == 'NovaSeqX' ):
        run_stats = get_run_info_novaseqX( instrument_model, application_version, tree, pipeline_type )
    elif( instrument_model == 'NovaSeq' ):
        run_stats = get_run_info_novaseq( instrument_model, application_version, tree, pipeline_type )
    elif( instrument_model == 'HiSeq3000' or instrument_model == 'HiSeq4000' ):
        run_stats = get_run_info_hiseq( instrument_model, application_version, tree )
    else:
        run_stats['instrument_type'] = 'unknown'

    return( run_stats )


def get_run_info_nextseq500( instrument_model, application_version, tree ):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        tree: xml tree
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    setup_node = tree.getroot().find("Setup")
    if setup_node is None:
        setup_node = tree.getroot()

    # Get required tree nodes.
    flowcell_node = tree.getroot().find("FlowCellRfidTag")

    # Now actually populate various stats
    run_stats['instrument_type'] = instrument_model
    run_stats['application_version'] = application_version
    run_stats['flow_cell_id'] = flowcell_node.find('SerialNumber').text
    run_stats['date'] = tree.getroot().find('RunStartDate').text
    run_stats['instrument'] = tree.getroot().find('InstrumentID').text
    run_stats['lanes'] = int(setup_node.find('NumLanes').text)
    run_stats['run_id'] = tree.getroot().find('RunID').text
    run_stats['r1_length'] = int(setup_node.find('Read1').text)
    run_stats['p7_index_length'] = int(setup_node.find('Index1Read').text)

    if( setup_node.find('Read1') != None ):
        run_stats['r2_length'] = int(setup_node.find('Read2').text)
        run_stats['p5_index_length'] = int(setup_node.find('Index2Read').text)
        run_stats['paired_end'] = True
    else:
        run_stats['paired_end'] = False

    run_stats['reverse_complement_i5'] = True

    return run_stats


def get_run_info_nextseq2000( instrument_model, application_version, tree ):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        tree: xml tree
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    setup_node = tree.getroot().find("Setup")
    if setup_node is None:
        setup_node = tree.getroot()

    # Get required tree nodes.
    completed_cycles = tree.getroot().find('CompletedCycles')

    # Now actually populate various stats
    run_stats['instrument_type'] = instrument_model
    run_stats['application_version'] = application_version
    run_stats['flow_cell_id'] = setup_node.find('FlowCellSerialNumber').text
    run_stats['date'] = setup_node.find('RunStartTime').text
    run_stats['instrument'] = setup_node.find('InstrumentSerialNumber').text
    run_stats['lanes'] = 1 # NextSeq 1000/2000 flowcells have one lane.

    output_dir = tree.getroot().find('OutputFolder').text
    run_stats['run_id'] = output_dir.split('/')[-1]

    run_stats['r1_length'] = int(completed_cycles.find('Read1').text)
    run_stats['p7_index_length'] = int(completed_cycles.find('Index1').text)

    if( completed_cycles.find('Read2') != None ):
        run_stats['r2_length'] = int(completed_cycles.find('Read2').text)
        run_stats['p5_index_length'] = int(completed_cycles.find('Index2').text)
        run_stats['paired_end'] = True
    else:
        run_stats['paired_end'] = False

    run_stats['reverse_complement_i5'] = True

    return run_stats


def get_run_info_miseq( instrument_model, application_version, tree ):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        tree: xml tree
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    setup_node = tree.getroot().find("Setup")
    if setup_node is None:
        setup_node = tree.getroot()

    # Get required tree nodes.
    flowcell_node = tree.getroot().find("FlowcellRFIDTag")
    reads_node = tree.getroot().find('Reads')

    # Now actually populate various stats
    run_stats['instrument_type'] = instrument_model
    run_stats['application_version'] = application_version
    run_stats['flow_cell_id'] = flowcell_node.find('SerialNumber').text
    run_stats['date'] = tree.getroot().find('RunStartDate').text
    run_stats['instrument'] = tree.getroot().find('ScannerID').text
    run_stats['lanes'] = int(setup_node.find('NumLanes').text)
    run_stats['run_id'] = tree.getroot().find('RunID').text

    read_len = []
    index_len = []
    for read_info in reads_node.findall('RunInfoRead'):
        attrib = read_info.attrib
        if( attrib['IsIndexedRead'] == 'Y' ):
          index_len.append( int( attrib['NumCycles'] ) )
        else:
          read_len.append( int( attrib['NumCycles'] ) )

    run_stats['r1_length'] = read_len[0]
    run_stats['p7_index_length'] = index_len[0]

    run_stats['paired_end'] = False
    if( len( read_len ) == 2 ):
      run_stats['r1_length'] = read_len[1]
      run_stats['p5_index_length'] = index_len[1]
      run_stats['paired_end'] = True

    run_stats['reverse_complement_i5'] = False

    return run_stats


def get_run_info_novaseq( instrument_model, application_version, tree, pipeline_type ):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        tree: xml tree
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    setup_node = tree.getroot().find("Setup")
    if setup_node is None:
        setup_node = tree.getroot()

    flowcell_node = tree.getroot().find("RfidsInfo")
    instrument_id_node = tree.getroot().find('InstrumentName')
    run_start_date_node = tree.getroot().find('RunStartDate')

    # Now actually populate various stats
    run_stats['instrument_type'] = instrument_model
    run_stats['application_version'] = application_version
    run_stats['flow_cell_id'] = flowcell_node.find('FlowCellSerialBarcode').text
    run_stats['date'] = run_start_date_node.text
    run_stats['instrument'] = instrument_id_node.text
    run_stats['flow_cell_mode'] = flowcell_node.find('FlowCellMode').text
    if( run_stats['flow_cell_mode'] in [ 'SP', 'S1', 'S2' ] ):
        run_stats['lanes'] = 2
    elif( run_stats['flow_cell_mode'] in [ 'S4' ] ):
        run_stats['lanes'] = 4
    else:
        raise ValueError( 'Unrecognized flow cell mode \'%s\'' % ( run_stats['flow_cell_mode'] ) )
    run_stats['run_id'] = tree.getroot().find('RunId').text

    # Read1 and Read2 may be absent
    run_stats['r1_length'] = int(setup_node.find('Read1NumberOfCycles').text)
    run_stats['p7_index_length'] = int(setup_node.find('IndexRead1NumberOfCycles').text)

    if( setup_node.find('Read2NumberOfCycles') != None ):
        run_stats['r2_length'] = int(setup_node.find('Read2NumberOfCycles').text)
        run_stats['p5_index_length'] = int(setup_node.find('IndexRead2NumberOfCycles').text)
        run_stats['paired_end'] = True
    else:
        run_stats['paired_end'] = False

    application = setup_node.find('Application').text
    application_version = setup_node.find('ApplicationVersion').text

    # Notes:
    #   o  NovaSeq application 1.7.0 can run reagent kit version 1.0 and 1.5
    #   o  the NWGC tells us:
    #      The NovaSeq v1.5 reagents are run on the NovaSeq that has an updated
    #      software which is version 1.7 that flips the i5 indices already on
    #      the sequencer when the data comes off. Typically, when the data
    #      comes off the sequencers, we need to flip both the i7 and i5 indices
    #      to the reverse complement in order to run fastqs or demux the data.
    #      With this being the case, only the i7 will need to be reverse
    #      complemented typically when data comes off the v1.5 version.
    #   o  however, not reverse complementing v1.5 fastqs in demultiplexing
    #      gives 'normal' looking sample-specific fastq files so I do not
    #      reverse complement here but allow for the possiblity in future
    #      reagent kits.
    #   o  The SBS consumable version differs between the two kits.  The line is
    #      <SbsConsumableVersion>1</SbsConsumableVersion>
    #        Key
    #        1= v1.0 SBS Reagents
    #        3= v1.5 SBS Reagents
    if( application_version == '1.7.0' or application_version == '1.7.5' ):
        sbs_consumable_version = flowcell_node.find('SbsConsumableVersion').text
        if( sbs_consumable_version == '1' ):
            run_stats['reverse_complement_i5'] = False
        elif( sbs_consumable_version == '3' ):
            if( pipeline_type == 'RNA-seq' ):
              run_stats['reverse_complement_i5'] = True
            elif( pipeline_type == 'ATAC-seq' ):
              run_stats['reverse_complement_i5'] = False
            else:
              raise ValueError('Unrecognized pipeline_type value \'%s\'' % ( pipeline_type ))
    elif( application_version == '1.8.0' or application_version == '1.8.1'):
        if( pipeline_type == 'RNA-seq' ):
            run_stats['reverse_complement_i5'] = True
        elif( pipeline_type == 'ATAC-seq' ):
            run_stats['reverse_complement_i5'] = False
        else:
            raise ValueError('Unrecognized pipeline_type value \'%s\'' % ( pipeline_type ))
    else:
        if( pipeline_type == 'RNA-seq' ):
            run_stats['reverse_complement_i5'] = True
        elif( pipeline_type == 'ATAC-seq' ):
            run_stats['reverse_complement_i5'] = False
        else:
            raise ValueError('Unrecognized pipeline_type value \'%s\'' % ( pipeline_type ))

    return run_stats


def get_run_info_novaseqX( instrument_model, application_version, tree, pipeline_type ):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        tree: xml tree
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    run_id_node = tree.getroot().find('RunId')
    instrument_id_node = tree.getroot().find('InstrumentSerialNumber')

    # Gather ConsumableInfo nodes.
    for consumable_info_outer in tree.getroot().iter('ConsumableInfo'):
        for consumable_info_inner in consumable_info_outer:
            for consumable_element in consumable_info_inner.iter():
                if(len(consumable_element) == 0):
                    continue
                if(consumable_element.find('Type').text == 'Buffer'):
                    consumable_element_buffer = consumable_element
                elif(consumable_element.find('Type').text == 'SampleTube'):
                    consumable_element_sample_tube = consumable_element
                elif(consumable_element.find('Type').text == 'Lyo'):
                    consumable_element_sample_lyo = consumable_element
                elif(consumable_element.find('Type').text == 'Reagent'):
                    consumable_element_sample_reagent = consumable_element
                elif(consumable_element.find('Type').text == 'FlowCell'):
                    consumable_element_sample_flowcell = consumable_element

    flowcell_node = consumable_element_sample_flowcell.find('SerialNumber').text

    # Now actually populate various stats
    run_stats['instrument_type'] = instrument_model
    run_stats['application_version'] = application_version
    run_stats['flow_cell_id'] = consumable_element_sample_flowcell.find('SerialNumber').text
    run_stats['date'] = run_id_node.text.split('_')[0]
    run_stats['instrument'] = instrument_id_node.text
    run_stats['flow_cell_mode'] = consumable_element_sample_flowcell.find('Mode').text
    run_stats['lanes'] = int(consumable_element_sample_tube.find('Name').text.split(' ')[0])

    run_stats['run_id'] = run_id_node.text

    run_stats['r1_length'] = None
    run_stats['r2_length'] = None
    run_stats['p7_index_length'] = None
    run_stats['p5_index_length'] = None
    run_stats['paired_end'] = False
    planned_reads_node = tree.getroot().find('PlannedReads')
    for read_node in planned_reads_node.iter('Read'):
        if(read_node.attrib['ReadName'] == 'Read1'):
            run_stats['r1_length'] = int(read_node.attrib['Cycles'])
        elif(read_node.attrib['ReadName'] == 'Read2'):
            run_stats['r2_length'] = int(read_node.attrib['Cycles'])
            run_stats['paired_end'] = True
        elif(read_node.attrib['ReadName'] == 'Index1'):
            run_stats['p7_index_length'] = int(read_node.attrib['Cycles'])
        elif(read_node.attrib['ReadName'] == 'Index2'):
            run_stats['p5_index_length'] = int(read_node.attrib['Cycles'])

    if( pipeline_type == 'RNA-seq' ):
        run_stats['reverse_complement_i5'] = True
    elif( pipeline_type == 'ATAC-seq' ):
        run_stats['reverse_complement_i5'] = False
    else:
        raise ValueError('Unrecognized pipeline_type value \'%s\'' % ( pipeline_type ))

    return run_stats


def get_run_info_hiseq( instrument_model, application_version, tree ):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        tree: xml tree
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    run_stats = {}

    setup_node = tree.getroot().find("Setup")
    if setup_node is None:
        setup_node = tree.getroot()

    run_id = setup_node.find('RunID').text

    # Now actually populate various stats
    run_stats['instrument_type'] = instrument_model
    run_stats['application_version'] = application_version
    run_stats['flow_cell_id'] = run_id.split('_')[3]
    run_stats['date'] = setup_node.find('RunStartDate').text
    run_stats['instrument'] = setup_node.find('ScannerID').text
    run_stats['lanes'] = 1
    run_stats['run_id'] = setup_node.find('RunID').text

    run_stats['r1_length'] = int(setup_node.find('Read1').text)
    run_stats['p7_index_length'] = int(setup_node.find('IndexRead1').text)

    if( setup_node.find('Read2') != None ):
        run_stats['r2_length'] = int(setup_node.find('Read2').text)
        run_stats['p5_index_length'] = int(setup_node.find('IndexRead2').text)
        run_stats['paired_end'] = True
    else:
        run_stats['paired_end'] = False

    if( instrument_model == 'HiSeq3000' ):
        run_stats['reverse_complement_i5'] = False
    elif( instrument_model == 'HiSeq4000' ):
        run_stats['reverse_complement_i5'] = True
    else:
        print( 'Error: unrecognized HiSeq model: \'%s\'' % ( instrument_model ) )
        sys.exit( -1 )

    return( run_stats )


# The following code is no longer used, I believe.
#
# def reverse_complement_i5(name):
#     """
#     Take a BCL directory or instrument type (NextSeq500, NextSeq2000, MiSeq,
#     NovaSeq, HiSeq4000, HiSeq3000, ...) and return whether or not i5 should
#     be reverse complemented. This assumes that NextSeq instruments and other
#     similar machines should be reverse complemented whereas MiSeq should not.
#     Args:
#         name (str): BCL directory or one of the instrument types as mentioned above
#
#     Returns:
#         bool: True if user would typically reverse complement i5 index and False otherwise.
#     """
#
#     if name in SEQUENCERS_P5_RC_MAP:
#         sequencer_type = name
#     elif os.path.exists(name):
#         sequencer_type = get_run_info(name)['instrument_type']
#
#         if sequencer_type not in SEQUENCERS_P5_RC_MAP:
#             raise ValueError('Sequencer type detected from BCL is %s, which is not in our known list of which sequencers require P5 reverse complementing or not.' % sequencer_type)
#     else:
#         raise ValueError('Invalid input, could not detect BCL or instrument ID.')
#
#     return SEQUENCERS_P5_RC_MAP[sequencer_type]
