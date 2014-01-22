from django.db import models

class RestrictionEnzyme(models.Model):
    name = models.CharField(max_length=15, unique=True)
    prototype = models.CharField(max_length=15)
    recognition_sequence = models.CharField(max_length=30)
    clean_recognition_sequence = models.CharField(max_length=30)
    suppliers = models.CharField(max_length=20)

    def __str__(self):
        return '{} {} {}'.format(self.name, self.recognition_sequence, self.suppliers)

