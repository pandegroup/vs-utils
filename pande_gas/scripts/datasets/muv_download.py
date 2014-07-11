#!/usr/bin/env python
"""
Download MUV datasets from PubChem.

The Java downloader from Rohrer and Baumann is broken, so we extract the
PubChem compound IDs (CIDs) from the jar and use the PubChem Power User
Gateway (PUG) to download compounds.

The jar is downloaded from the website of the MUV authors:
http://www.pharmchem.tu-bs.de/lehre/baumann/MUV.html

See Rohrer and Baumann, Maximum Unbiased Validation (MUV) Data Sets for
Virtual Screening Based on PubChem Bioactivity Data. J. Chem. Inf. Model.,
2009, 49 (2), pp 169-184.
"""

__author__ = "Steven Kearnes"
__copyright__ = "Copyright 2014, Stanford University"
__license__ = "BSD 3-clause"

import os
import re
import sys
import tempfile
import time
import urllib
import urllib2
import zipfile


def main():
    muv_url = 'http://www.pharmchem.tu-bs.de/lehre/baumann/MUV_Downloader.jar'
    jar, _ = urllib.urlretrieve(muv_url)
    with zipfile.ZipFile(jar) as f:
        path = tempfile.mkdtemp()
        f.extract('main/main.jar', path)
    with zipfile.ZipFile(path + '/main/main.jar') as f:
        for m in f.filelist:
            basename = os.path.basename(m.filename)
            if basename.startswith('aid'):
                output = os.path.splitext(basename)[0] + '.sdf.gz'
                print '{} -> {}'.format(basename, output)
                with f.open(m) as g:
                    cids = [line.strip() for line in g]
                    pc_url = cid_query(cids)
                    urllib.urlretrieve(pc_url, output)
                    print 'done.\n'


def cid_query(cids):
    """
    Download a batch of CIDs from PubChem using PUG.

    Parameters
    ----------
    cids : list
        PubChem compound IDs (CIDs).
    """
    pug = 'https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi'
    query_template = """
<PCT-Data>
 <PCT-Data_input>
  <PCT-InputData>
   <PCT-InputData_download>
    <PCT-Download>
     <PCT-Download_uids>
      <PCT-QueryUids>
       <PCT-QueryUids_ids>
        <PCT-ID-List>
         <PCT-ID-List_db>pccompound</PCT-ID-List_db>
         <PCT-ID-List_uids>
          %(uids)s
         </PCT-ID-List_uids>
        </PCT-ID-List>
       </PCT-QueryUids_ids>
      </PCT-QueryUids>
     </PCT-Download_uids>
     <PCT-Download_format value="sdf"/>
     <PCT-Download_compression value="gzip"/>
    </PCT-Download>
   </PCT-InputData_download>
  </PCT-InputData>
 </PCT-Data_input>
</PCT-Data>
"""
    status_template = """
<PCT-Data>
 <PCT-Data_input>
  <PCT-InputData>
   <PCT-InputData_request>
    <PCT-Request>
     <PCT-Request_reqid>%(wait_id)s</PCT-Request_reqid>
     <PCT-Request_type value="status"/>
    </PCT-Request>
   </PCT-InputData_request>
  </PCT-InputData>
 </PCT-Data_input>
</PCT-Data>
"""
    xml_cids = ''
    for cid in cids:
        xml_cids += '<PCT-ID-List_uids_E>{}</PCT-ID-List_uids_E>\n'.format(cid)
    query = query_template % {'uids': xml_cids}
    url = None
    while True:
        q = urllib2.urlopen(pug, query)
        response = q.read()
        url_re = re.search(
            '<PCT-Download-URL_url>\s*(.*?)\s*</PCT-Download-URL_url>',
            response)
        if url_re is not None:
            url = url_re.groups()[0]
            break
        sys.stdout.write('.')
        sys.stdout.flush()
        wait_re = re.search(
            '<PCT-Waiting_reqid>\s*(.*?)\s*</PCT-Waiting_reqid>', response)
        wait_id = wait_re.groups()[0]
        query = status_template % {'wait_id': wait_id}
        time.sleep(10)
    return url

if __name__ == '__main__':
    main()
