from django.db import models


class PrototypeEnzyme(models.Model):
    name = models.CharField(max_length=15)
    clean_recognition_sequence = models.CharField(max_length=30)

    def __str__(self):
        res = '\n'
        for re in self.restrictionenzyme_set.all():
            res += '\t{} {} {}\n'.format(re.name, re.recognition_sequence, re.suppliers)
        return '{} {} {}'.format(self.name, self.clean_recognition_sequence,res)

class RestrictionEnzyme(models.Model):
    prototype = models.ForeignKey(PrototypeEnzyme)
    name = models.CharField(max_length=15)
    recognition_sequence = models.CharField(max_length=30)
    suppliers = models.CharField(max_length=20)

    def __str__(self):
        return '{} {} {}'.format(self.name, self.recognition_sequence, self.suppliers)



