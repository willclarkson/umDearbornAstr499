# Test Routines

# LDAC routines

# test = catalog('ldac test')
# read = test.read_fits('testdata/test.ldac', maxflag=0)
# print '%d sources and %d columns read from ldac file' % \
#     (read[0], read[1])

# print test.write_ascii('testdata/test.ldac.dat'), \
#     'sources sucessfully written to ascii file'

# print test.write_fits('testdata/test_written.ldac'), \
#     'sources successfully written to fits file'

# print test.write_database('testdata/test.ldac.db')


# catalog download and manipulation

# cat1 = catalog('APASS9')
# print(cat1.download_catalog(80, 0, 0.5, 10000), 'sources grabbed from', cat1.catalogname)
# print(cat1.fields)
# print(cat1[0])

# cat2 = catalog('2MASS')
# print cat2.download_catalog(80, 0, 0.5, 10000), 'sources grabbed from', cat2.catalogname
# print cat2[305]
# print cat2.fields

# cat3 = catalog('SDSS-R9')
# print(cat3.download_catalog(329.50922672, -2.33703204, 0.001, 10000), 'sources grabbed from', cat3.catalogname)
# print(cat3.data['imag'])


# cat4 = catalog('USNO-B1')
# print cat4.download_catalog(80, 0, 0.5, 10000), 'sources grabbed from', cat4.catalogname
# print cat4[305]
# print cat4.fields


# print cat2.shape, cat2.history
# print cat2.transform_filters('V'), 'sources transformed to V'      #2MASS to BVRI
# print cat2.shape, cat2.history

# print cat2.shape, cat2.history
# print cat2.transform_filters('K_UKIRT'), 'sources transformed to K'      #2MASS to UKIRT
# print cat2.shape, cat2.history

# print cat3.shape, cat3.history
# print cat3.transform_filters('Z_UKIRT'), 'sources transformed to Z'      #SDSS to UKIRT
# print cat3.shape, cat3.history

# print cat1.shape, cat1.history
# print cat1.transform_filters('I'), 'sources transformed to I'      #SDSS to BVRI
# print cat1.shape, cat1.history
# for i in cat1:
#    print i['_Imag'], i['rmag'], i['imag']


# print test.write_database('test.db'), 'sources written to database file'

# cat4 = catalog('')
# print cat4.read_database('test.db'), 'sources read from database file'


# test Gaia
cat = catalog('GAIA')
print(cat.download_catalog(294.99525, 0.065194, 0.5, 10000),
      'sources grabbed from ', cat.catalogname)

print(cat.data)


# compare GAIA against SDSS
cat = catalog('GAIA')
print(cat.download_catalog(329.50922672, -2.33703204, 0.1, 10000),
      'sources grabbed from', cat.catalogname)
print(cat.transform_filters('I'), 'stars transformed')

refcat = catalog('SDSS-R9')
print(refcat.download_catalog(329.50922672, -2.33703204, 0.1, 10000),
      'sources grabbed from', refcat.catalogname)
refcat.transform_filters('I')

refmag, gaiamag = refcat.match_with(cat, extract_this_catalog=['_Imag'],
                                    extract_other_catalog=['_Imag'])

print(np.median(np.array(refmag)-np.array(gaiamag)),
      np.std(np.array(refmag)-np.array(gaiamag)))
